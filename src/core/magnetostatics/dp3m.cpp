/*
 * Copyright (C) 2010-2024 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 *  P3M algorithm for long-range magnetic dipole-dipole interaction.
 *
 *  By default the magnetic epsilon is metallic = 0.
 */

#include "config/config.hpp"

#ifdef DP3M

#include "magnetostatics/dp3m.hpp"
#include "magnetostatics/dp3m.impl.hpp"

#include "fft/fft.hpp"
#include "p3m/TuningAlgorithm.hpp"
#include "p3m/TuningLogger.hpp"
#include "p3m/common.hpp"
#include "p3m/influence_function_dipolar.hpp"
#include "p3m/interpolation.hpp"
#include "p3m/math.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrators/Propagation.hpp"
#include "npt.hpp"
#include "system/System.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iterator>
#include <numbers>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#ifdef FFTW3_H
#error "The FFTW3 library shouldn't be visible in this translation unit"
#endif

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::count_magnetic_particles() {
  int local_n = 0;
  double local_mu2 = 0.;

  for (auto const &p : get_system().cell_structure->local_particles()) {
    if (p.dipm() != 0.) {
      local_mu2 += p.calc_dip().norm2();
      local_n++;
    }
  }

  boost::mpi::all_reduce(comm_cart, local_mu2, dp3m.sum_mu2, std::plus<>());
  boost::mpi::all_reduce(comm_cart, local_n, dp3m.sum_dip_part, std::plus<>());
}

static double dp3m_k_space_error(double box_size, int mesh, int cao,
                                 int n_c_part, double sum_q2, double alpha_L);

static double dp3m_real_space_error(double box_size, double r_cut_iL,
                                    int n_c_part, double sum_q2,
                                    double alpha_L);

/** Compute the value of alpha through a bisection method.
 *  Based on eq. (33) @cite wang01a.
 */
double dp3m_rtbisection(double box_size, double r_cut_iL, int n_c_part,
                        double sum_q2, double x1, double x2, double xacc,
                        double tuned_accuracy);

template <typename FloatType, Arch Architecture>
double
DipolarP3MImpl<FloatType, Architecture>::calc_average_self_energy_k_space()
    const {
  auto const &box_geo = *get_system().box_geo;
  auto const node_phi = grid_influence_function_self_energy(
      dp3m.params, dp3m.mesh.start, dp3m.mesh.stop, dp3m.g_energy);

  double phi = 0.;
  boost::mpi::reduce(comm_cart, node_phi, phi, std::plus<>(), 0);
  phi /= 3. * box_geo.length()[0] * Utils::int_pow<3>(dp3m.params.mesh[0]);
  return phi * std::numbers::pi;
}

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::init_cpu_kernels() {
  assert(dp3m.params.mesh >= Utils::Vector3i::broadcast(1));
  assert(dp3m.params.cao >= 1 and dp3m.params.cao <= 7);
  assert(dp3m.params.alpha > 0.);

  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  auto const verlet_skin = system.cell_structure->get_verlet_skin();

  dp3m.params.cao3 = Utils::int_pow<3>(dp3m.params.cao);
  dp3m.params.recalc_a_ai_cao_cut(box_geo.length());

  assert(dp3m.fft);
  dp3m.local_mesh.calc_local_ca_mesh(dp3m.params, local_geo, verlet_skin, 0.);
  dp3m.fft_buffers->init_halo();
  dp3m.fft->init(dp3m.params);
  dp3m.mesh.ks_pnum = dp3m.fft->get_ks_pnum();
  dp3m.fft_buffers->init_meshes(dp3m.fft->get_ca_mesh_size());
  dp3m.update_mesh_views();
  dp3m.calc_differential_operator();

  /* fix box length dependent constants */
  scaleby_box_l();

  count_magnetic_particles();
}

namespace {
template <int cao> struct AssignDipole {
  void operator()(auto &dp3m, Utils::Vector3d const &real_pos,
                  Utils::Vector3d const &dip) const {
    using value_type =
        typename std::remove_reference_t<decltype(dp3m)>::value_type;
    auto const weights = p3m_calculate_interpolation_weights<cao>(
        real_pos, dp3m.params.ai, dp3m.local_mesh);
    p3m_interpolate<cao>(
        dp3m.local_mesh, weights, [&dip, &dp3m](int ind, double w) {
          dp3m.mesh.rs_fields[0u][ind] += value_type(w * dip[0u]);
          dp3m.mesh.rs_fields[1u][ind] += value_type(w * dip[1u]);
          dp3m.mesh.rs_fields[2u][ind] += value_type(w * dip[2u]);
        });

    dp3m.inter_weights.template store<cao>(weights);
  }
};
} // namespace

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::dipole_assign(
    ParticleRange const &particles) {
  dp3m.inter_weights.reset(dp3m.params.cao);

  /* prepare local FFT mesh */
  for (auto &rs_mesh_field : dp3m.mesh.rs_fields)
    for (int j = 0; j < dp3m.local_mesh.size; j++)
      rs_mesh_field[j] = 0.;

  for (auto const &p : particles) {
    if (p.dipm() != 0.) {
      Utils::integral_parameter<int, AssignDipole, 1, 7>(dp3m.params.cao, dp3m,
                                                         p.pos(), p.calc_dip());
    }
  }
}

namespace {
template <int cao> struct AssignTorques {
  void operator()(auto &dp3m, double prefac, int d_rs,
                  ParticleRange const &particles) const {

    /* magnetic particle index */
    auto p_index = std::size_t{0ul};

    for (auto &p : particles) {
      if (p.dipm() != 0.) {
        auto const weights = dp3m.inter_weights.template load<cao>(p_index);

        Utils::Vector3d E{};
        p3m_interpolate(dp3m.local_mesh, weights,
                        [&E, &dp3m, d_rs](int ind, double w) {
                          E[d_rs] += w * double(dp3m.mesh.rs_scalar[ind]);
                        });

        p.torque() -= vector_product(p.calc_dip(), prefac * E);
        ++p_index;
      }
    }
  }
};

template <int cao> struct AssignForces {
  void operator()(auto &dp3m, double prefac, int d_rs,
                  ParticleRange const &particles) const {

    /* magnetic particle index */
    auto p_index = std::size_t{0ul};

    for (auto &p : particles) {
      if (p.dipm() != 0.) {
        auto const weights = dp3m.inter_weights.template load<cao>(p_index);

        Utils::Vector3d E{};
        p3m_interpolate(dp3m.local_mesh, weights,
                        [&E, &dp3m](int ind, double w) {
                          E[0u] += w * double(dp3m.mesh.rs_fields[0u][ind]);
                          E[1u] += w * double(dp3m.mesh.rs_fields[1u][ind]);
                          E[2u] += w * double(dp3m.mesh.rs_fields[2u][ind]);
                        });

        p.force()[d_rs] += p.calc_dip() * prefac * E;
        ++p_index;
      }
    }
  }
};
} // namespace

template <typename FloatType, Arch Architecture>
double DipolarP3MImpl<FloatType, Architecture>::long_range_kernel(
    bool force_flag, bool energy_flag, ParticleRange const &particles) {
  /* k-space energy */
  double energy = 0.;
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const dipole_prefac = prefactor / Utils::int_pow<3>(dp3m.params.mesh[0]);
#ifdef NPT
  auto const npt_flag =
      force_flag and (system.propagation->integ_switch == INTEG_METHOD_NPT_ISO);
#else
  auto constexpr npt_flag = false;
#endif

  if (dp3m.sum_mu2 > 0.) {
    dipole_assign(particles);
    dp3m.fft_buffers->perform_vector_halo_gather();
    for (auto &rs_mesh : dp3m.fft_buffers->get_vector_mesh()) {
      dp3m.fft->forward_fft(rs_mesh);
    }
    dp3m.update_mesh_views();
  }

  /* === k-space energy calculation  === */
  if (energy_flag or npt_flag) {
    /*********************
       Dipolar energy
    **********************/
    if (dp3m.sum_mu2 > 0.) {
      /* i*k differentiation for dipolar gradients:
       * |(\Fourier{\vect{mu}}(k)\cdot \vect{k})|^2 */

      auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
      auto const &offset = dp3m.mesh.start;
      auto const &d_op = dp3m.d_op[0u];
      auto const &mesh_dip = dp3m.mesh.rs_fields;
      auto const [KX, KY, KZ] = dp3m.fft->get_permutations();
      auto indices = Utils::Vector3i{};
      auto index = std::size_t(0u);
      auto it_energy = dp3m.g_energy.begin();
      auto node_energy = 0.;
      for_each_3d(mesh_start, dp3m.mesh.size, indices, [&]() {
        auto const shift = indices + offset;
        // Re(mu)*k
        auto const re = mesh_dip[0u][index] * FloatType(d_op[shift[KX]]) +
                        mesh_dip[1u][index] * FloatType(d_op[shift[KY]]) +
                        mesh_dip[2u][index] * FloatType(d_op[shift[KZ]]);
        ++index;
        // Im(mu)*k
        auto const im = mesh_dip[0u][index] * FloatType(d_op[shift[KX]]) +
                        mesh_dip[1u][index] * FloatType(d_op[shift[KY]]) +
                        mesh_dip[2u][index] * FloatType(d_op[shift[KZ]]);
        ++index;
        node_energy += *it_energy * (Utils::sqr(re) + Utils::sqr(im));
        std::advance(it_energy, 1);
      });

      node_energy *= dipole_prefac * std::numbers::pi * box_geo.length_inv()[0];
      boost::mpi::reduce(comm_cart, node_energy, energy, std::plus<>(), 0);

      if (dp3m.energy_correction == 0.)
        calc_energy_correction();

      if (this_node == 0) {
        /* self energy correction */
        energy -= prefactor * dp3m.sum_mu2 * std::numbers::inv_sqrtpi *
                  (2. / 3.) * Utils::int_pow<3>(dp3m.params.alpha);

        /* dipolar energy correction due to systematic Madelung-self effects */
        energy += prefactor * dp3m.energy_correction / box_geo.volume();
      }
    }
  } // if (energy_flag)

  /* === k-space force calculation  === */
  if (force_flag) {
    /****************************
     * DIPOLAR TORQUES (k-space)
     ****************************/
    if (dp3m.sum_mu2 > 0.) {
      auto const wavenumber = 2. * std::numbers::pi * box_geo.length_inv()[0u];
      auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
      auto const &mesh_stop = dp3m.mesh.size;
      auto const offset = dp3m.mesh.start;
      auto const &d_op = dp3m.d_op[0u];
      auto const [KX, KY, KZ] = dp3m.fft->get_permutations();
      auto &mesh_dip = dp3m.mesh.rs_fields;
      auto indices = Utils::Vector3i{};
      auto index = std::size_t(0u);
      dp3m.ks_scalar.resize(dp3m.mesh.rs_scalar.size());

      /* fill in ks_scalar array for torque calculation */
      auto it_energy = dp3m.g_energy.begin();
      index = 0u;
      for_each_3d(mesh_start, mesh_stop, indices, [&]() {
        auto const shift = indices + offset;
        // Re(mu)*k
        auto const re = mesh_dip[0u][index] * FloatType(d_op[shift[KX]]) +
                        mesh_dip[1u][index] * FloatType(d_op[shift[KY]]) +
                        mesh_dip[2u][index] * FloatType(d_op[shift[KZ]]);
        dp3m.ks_scalar[index] = *it_energy * re;
        ++index;
        // Im(mu)*k
        auto const im = mesh_dip[0u][index] * FloatType(d_op[shift[KX]]) +
                        mesh_dip[1u][index] * FloatType(d_op[shift[KY]]) +
                        mesh_dip[2u][index] * FloatType(d_op[shift[KZ]]);
        dp3m.ks_scalar[index] = *it_energy * im;
        ++index;
        std::advance(it_energy, 1);
      });

      /* Force component loop */
      for (int d = 0; d < 3; d++) {
        index = 0u;
        for_each_3d(mesh_start, mesh_stop, indices, [&]() {
          auto const d_op_val = FloatType(d_op[indices[d] + offset[d]]);
          dp3m.mesh.rs_scalar[index] = d_op_val * dp3m.ks_scalar[index];
          ++index;
          dp3m.mesh.rs_scalar[index] = d_op_val * dp3m.ks_scalar[index];
          ++index;
        });
        dp3m.fft->backward_fft(dp3m.fft_buffers->get_scalar_mesh());
        dp3m.fft_buffers->perform_scalar_halo_spread();
        /* Assign force component from mesh to particle */
        auto const d_rs = (d + dp3m.mesh.ks_pnum) % 3;
        Utils::integral_parameter<int, AssignTorques, 1, 7>(
            dp3m.params.cao, dp3m, dipole_prefac * wavenumber, d_rs, particles);
      }

      /***************************
         DIPOLAR FORCES (k-space)
      ****************************/
      // Compute forces after torques because the algorithm below overwrites the
      // grids dp3m.mesh.rs_fields !
      // Note: I'll do here 9 inverse FFTs. By symmetry, we can reduce this
      // number to 6 !
      /* fill in ks_scalar array for force calculation */
      auto it_force = dp3m.g_force.begin();
      index = 0u;
      for_each_3d(mesh_start, mesh_stop, indices, [&]() {
        auto const shift = indices + offset;
        // Re(mu)*k
        auto const re = mesh_dip[0u][index] * FloatType(d_op[shift[KX]]) +
                        mesh_dip[1u][index] * FloatType(d_op[shift[KY]]) +
                        mesh_dip[2u][index] * FloatType(d_op[shift[KZ]]);
        ++index;
        // Im(mu)*k
        auto const im = mesh_dip[0u][index] * FloatType(d_op[shift[KX]]) +
                        mesh_dip[1u][index] * FloatType(d_op[shift[KY]]) +
                        mesh_dip[2u][index] * FloatType(d_op[shift[KZ]]);
        ++index;
        dp3m.ks_scalar[index - 2] = *it_force * im;
        dp3m.ks_scalar[index - 1] = *it_force * (-re);
        std::advance(it_force, 1);
      });

      /* Force component loop */
      for (int d = 0; d < 3; d++) {
        index = 0u;
        for_each_3d(mesh_start, mesh_stop, indices, [&]() {
          auto const d_op_val = FloatType(d_op[indices[d] + offset[d]]);
          auto const shift = indices + offset;
          auto const f1 = d_op_val * dp3m.ks_scalar[index];
          mesh_dip[0u][index] = FloatType(d_op[shift[KX]]) * f1;
          mesh_dip[1u][index] = FloatType(d_op[shift[KY]]) * f1;
          mesh_dip[2u][index] = FloatType(d_op[shift[KZ]]) * f1;
          ++index;
          auto const f2 = d_op_val * dp3m.ks_scalar[index];
          mesh_dip[0u][index] = FloatType(d_op[shift[KX]]) * f2;
          mesh_dip[1u][index] = FloatType(d_op[shift[KY]]) * f2;
          mesh_dip[2u][index] = FloatType(d_op[shift[KZ]]) * f2;
          ++index;
        });
        for (auto &rs_mesh : dp3m.fft_buffers->get_vector_mesh()) {
          dp3m.fft->backward_fft(rs_mesh);
        }
        dp3m.fft_buffers->perform_vector_halo_spread();
        /* Assign force component from mesh to particle */
        auto const d_rs = (d + dp3m.mesh.ks_pnum) % 3;
        Utils::integral_parameter<int, AssignForces, 1, 7>(
            dp3m.params.cao, dp3m, dipole_prefac * Utils::sqr(wavenumber), d_rs,
            particles);
      }
    } /* if (dp3m.sum_mu2 > 0) */
  } /* if (force_flag) */

  if (dp3m.params.epsilon != P3M_EPSILON_METALLIC) {
    auto const surface_term =
        calc_surface_term(force_flag, energy_flag or npt_flag, particles);
    if (this_node == 0) {
      energy += surface_term;
    }
  }
  if (npt_flag) {
    npt_add_virial_contribution(energy);
    fprintf(stderr, "dipolar_P3M at this moment is added to p_vir[0]\n");
  }
  if (not energy_flag) {
    energy = 0.;
  }

  return energy;
}

template <typename FloatType, Arch Architecture>
double DipolarP3MImpl<FloatType, Architecture>::calc_surface_term(
    bool force_flag, bool energy_flag, ParticleRange const &particles) {
  auto const &box_geo = *get_system().box_geo;
  auto const pref = prefactor * 4. * std::numbers::pi / box_geo.volume() /
                    (2. * dp3m.params.epsilon + 1.);
  auto const n_local_part = particles.size();

  // We put all the dipolar momenta in a the arrays mx,my,mz according to the
  // id-number of the particles
  std::vector<double> mx(n_local_part);
  std::vector<double> my(n_local_part);
  std::vector<double> mz(n_local_part);

  std::size_t ip = 0u;
  for (auto const &p : particles) {
    auto const dip = p.calc_dip();
    mx[ip] = dip[0u];
    my[ip] = dip[1u];
    mz[ip] = dip[2u];
    ip++;
  }

  // we will need the sum of all dipolar momenta vectors
  auto local_dip = Utils::Vector3d{};
  for (std::size_t i = 0u; i < n_local_part; i++) {
    local_dip[0u] += mx[i];
    local_dip[1u] += my[i];
    local_dip[2u] += mz[i];
  }
  auto const box_dip =
      boost::mpi::all_reduce(comm_cart, local_dip, std::plus<>());

  double energy = 0.;
  if (energy_flag) {
    double sum_e = 0.;
    for (std::size_t i = 0u; i < n_local_part; i++) {
      sum_e += mx[i] * box_dip[0] + my[i] * box_dip[1] + mz[i] * box_dip[2];
    }
    energy =
        0.5 * pref * boost::mpi::all_reduce(comm_cart, sum_e, std::plus<>());
  }

  if (force_flag) {

    std::vector<double> sumix(n_local_part);
    std::vector<double> sumiy(n_local_part);
    std::vector<double> sumiz(n_local_part);

    for (std::size_t i = 0u; i < n_local_part; i++) {
      sumix[i] = my[i] * box_dip[2u] - mz[i] * box_dip[1u];
      sumiy[i] = mz[i] * box_dip[0u] - mx[i] * box_dip[2u];
      sumiz[i] = mx[i] * box_dip[1u] - my[i] * box_dip[0u];
    }

    ip = 0u;
    for (auto &p : particles) {
      auto &torque = p.torque();
      torque[0u] -= pref * sumix[ip];
      torque[1u] -= pref * sumiy[ip];
      torque[2u] -= pref * sumiz[ip];
      ip++;
    }
  }

  return energy;
}

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::calc_influence_function_force() {
  dp3m.g_force = grid_influence_function<FloatType, 3>(
      dp3m.params, dp3m.mesh.start, dp3m.mesh.stop,
      get_system().box_geo->length_inv());
}

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::calc_influence_function_energy() {
  dp3m.g_energy = grid_influence_function<FloatType, 2>(
      dp3m.params, dp3m.mesh.start, dp3m.mesh.stop,
      get_system().box_geo->length_inv());
}

template <typename FloatType, Arch Architecture>
class DipolarTuningAlgorithm : public TuningAlgorithm {
  p3m_data_struct_dipoles<FloatType> &dp3m;
  int m_mesh_max = -1, m_mesh_min = -1;

public:
  DipolarTuningAlgorithm(System::System &system, decltype(dp3m) &input_dp3m,
                         double prefactor, int timings)
      : TuningAlgorithm(system, prefactor, timings), dp3m{input_dp3m} {}

  P3MParameters &get_params() override { return dp3m.params; }

  void on_solver_change() const override { m_system.on_dipoles_change(); }

  std::optional<std::string>
  layer_correction_veto_r_cut(double) const override {
    return {};
  }

  void setup_logger(bool verbose) override {
    auto const &box_geo = *m_system.box_geo;
    m_logger = std::make_unique<TuningLogger>(
        verbose and this_node == 0, "DipolarP3M", TuningLogger::Mode::Dipolar);
    m_logger->tuning_goals(dp3m.params.accuracy, m_prefactor,
                           box_geo.length()[0], dp3m.sum_dip_part,
                           dp3m.sum_mu2);
    m_logger->log_tuning_start();
  }

  std::tuple<double, double, double, double>
  calculate_accuracy(Utils::Vector3i const &mesh, int cao,
                     double r_cut_iL) const override {

    double alpha_L, rs_err, ks_err;
    auto const &box_geo = *m_system.box_geo;

    /* calc maximal real space error for setting */
    rs_err = dp3m_real_space_error(box_geo.length()[0], r_cut_iL,
                                   dp3m.sum_dip_part, dp3m.sum_mu2, 0.001);
    // alpha cannot be zero for dipoles because real-space formula breaks down

    if (std::numbers::sqrt2 * rs_err > dp3m.params.accuracy) {
      /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
      alpha_L = dp3m_rtbisection(
          box_geo.length()[0], r_cut_iL, dp3m.sum_dip_part, dp3m.sum_mu2,
          0.0001 * box_geo.length()[0], 5. * box_geo.length()[0], 0.0001,
          dp3m.params.accuracy);
    } else {
      /* even alpha=0 is ok, however, we cannot choose it since it kills the
         k-space error formula.
         Anyways, this very likely NOT the optimal solution */
      alpha_L = 0.1;
    }

    /* calculate real-space and k-space error for this alpha_L */
    rs_err = dp3m_real_space_error(box_geo.length()[0], r_cut_iL,
                                   dp3m.sum_dip_part, dp3m.sum_mu2, alpha_L);
    ks_err = dp3m_k_space_error(box_geo.length()[0], mesh[0], cao,
                                dp3m.sum_dip_part, dp3m.sum_mu2, alpha_L);

    return {Utils::Vector2d{rs_err, ks_err}.norm(), rs_err, ks_err, alpha_L};
  }

  void determine_mesh_limits() override {
    if (dp3m.params.mesh[0] == -1) {
      /* simple heuristic to limit the tried meshes if the accuracy cannot
         be obtained with smaller meshes, but normally not all these
         meshes have to be tested */
      auto const expo = std::log(std::cbrt(dp3m.sum_dip_part)) / std::log(2.);
      /* Medium-educated guess for the minimal mesh */
      m_mesh_min = static_cast<int>(std::round(std::pow(2., std::floor(expo))));
      /* avoid using more than 1 GB of FFT arrays */
      m_mesh_max = 128;
    } else {
      m_mesh_min = m_mesh_max = dp3m.params.mesh[0];
      m_logger->report_fixed_mesh(dp3m.params.mesh);
    }
  }

  TuningAlgorithm::Parameters get_time() override {
    auto tuned_params = TuningAlgorithm::Parameters{};
    auto time_best = time_sentinel;
    for (auto tmp_mesh = m_mesh_min; tmp_mesh <= m_mesh_max; tmp_mesh += 2) {
      auto trial_params = TuningAlgorithm::Parameters{};
      trial_params.mesh = Utils::Vector3i::broadcast(tmp_mesh);
      trial_params.cao = cao_best;

      auto const trial_time =
          get_m_time(trial_params.mesh, trial_params.cao, trial_params.r_cut_iL,
                     trial_params.alpha_L, trial_params.accuracy);

      /* this mesh does not work at all */
      if (trial_time < 0.)
        continue;

      /* the optimum r_cut for this mesh is the upper limit for higher meshes,
         everything else is slower */
      m_r_cut_iL_max = trial_params.r_cut_iL;

      if (trial_time < time_best) {
        /* new optimum */
        reset_n_trials();
        tuned_params = trial_params;
        time_best = tuned_params.time = trial_time;
      } else if (trial_time > time_best + time_granularity or
                 get_n_trials() > max_n_consecutive_trials) {
        /* no hope of further optimisation */
        break;
      }
    }
    return tuned_params;
  }
};

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::tune() {
  auto &system = get_system();
  auto const &box_geo = *system.box_geo;
  if (dp3m.params.alpha_L == 0. and dp3m.params.alpha != 0.) {
    dp3m.params.alpha_L = dp3m.params.alpha * box_geo.length()[0];
  }
  if (dp3m.params.r_cut_iL == 0. and dp3m.params.r_cut != 0.) {
    dp3m.params.r_cut_iL = dp3m.params.r_cut * box_geo.length_inv()[0];
  }
  if (not is_tuned()) {
    count_magnetic_particles();
    if (dp3m.sum_dip_part == 0) {
      throw std::runtime_error(
          "DipolarP3M: no dipolar particles in the system");
    }
    try {
      DipolarTuningAlgorithm<FloatType, Architecture> parameters(
          system, dp3m, prefactor, tune_timings);
      parameters.setup_logger(tune_verbose);
      // parameter ranges
      parameters.determine_mesh_limits();
      parameters.determine_r_cut_limits();
      parameters.determine_cao_limits(3);
      // run tuning algorithm
      parameters.tune();
      m_is_tuned = true;
      system.on_dipoles_change();
    } catch (...) {
      dp3m.params.tuning = false;
      throw;
    }
  }
  init();
}

/** Tuning dipolar-P3M */
static auto dp3m_tune_aliasing_sums(Utils::Vector3i const &shift, int mesh,
                                    double mesh_i, int cao, double alpha_L_i) {

  auto constexpr mesh_start = Utils::Vector3i::broadcast(-P3M_BRILLOUIN);
  auto constexpr mesh_stop = Utils::Vector3i::broadcast(P3M_BRILLOUIN + 1);
  auto const factor1 = Utils::sqr(std::numbers::pi * alpha_L_i);
  auto alias1 = 0.;
  auto alias2 = 0.;

  Utils::Vector3i indices{};
  Utils::Vector3i nm{};
  Utils::Vector3d fnm{};
  for_each_3d(
      mesh_start, mesh_stop, indices,
      [&]() {
        auto const norm_sq = nm.norm2();
        auto const ex = std::exp(-factor1 * norm_sq);
        auto const U2 = std::pow(Utils::product(fnm), 2 * cao);
        alias1 += Utils::sqr(ex) * norm_sq;
        alias2 += U2 * ex * std::pow(shift * nm, 3) / norm_sq;
      },
      [&](unsigned dim, int n) {
        nm[dim] = shift[dim] + n * mesh;
        fnm[dim] = math::sinc(nm[dim] * mesh_i);
      });

  return std::make_pair(alias1, alias2);
}

/** Calculate the k-space error of dipolar-P3M */
static double dp3m_k_space_error(double box_size, int mesh, int cao,
                                 int n_c_part, double sum_q2, double alpha_L) {

  auto const cotangent_sum = math::get_analytic_cotangent_sum_kernel(cao);
  auto const mesh_i = 1. / static_cast<double>(mesh);
  auto const alpha_L_i = 1. / alpha_L;
  auto const mesh_stop = Utils::Vector3i::broadcast(mesh / 2);
  auto const mesh_start = -mesh_stop;
  auto indices = Utils::Vector3i{};
  auto values = Utils::Vector3d{};
  auto he_q = 0.;

  for_each_3d(
      mesh_start, mesh_stop, indices,
      [&]() {
        if ((indices[0] != 0) or (indices[1] != 0) or (indices[2] != 0)) {
          auto const n2 = indices.norm2();
          auto const cs = Utils::product(values);
          auto const [alias1, alias2] =
              dp3m_tune_aliasing_sums(indices, mesh, mesh_i, cao, alpha_L_i);
          auto const d =
              alias1 - Utils::sqr(alias2 / cs) /
                           Utils::int_pow<3>(static_cast<double>(n2));
          /* at high precision, d can become negative due to extinction;
             also, don't take values that have no significant digits left*/
          if (d > 0. and std::fabs(d / alias1) > ROUND_ERROR_PREC)
            he_q += d;
        }
      },
      [&values, &mesh_i, cotangent_sum](unsigned dim, int n) {
        values[dim] = cotangent_sum(n, mesh_i);
      });

  return 8. * Utils::sqr(std::numbers::pi) / 3. * sum_q2 *
         sqrt(he_q / n_c_part) / Utils::int_pow<4>(box_size);
}

/** Calculate the value of the errors for the REAL part of the force in terms
 *  of the splitting parameter alpha of Ewald. Based on eq. (33) @cite wang01a.
 *
 *  Please note that in this more refined approach we don't use
 *  eq. (37), but eq. (33) which maintains all the powers in alpha.
 */
static double dp3m_real_space_error(double box_size, double r_cut_iL,
                                    int n_c_part, double sum_q2,
                                    double alpha_L) {
  double d_error_f, d_cc, d_dc, d_con;

  auto const d_rcut = r_cut_iL * box_size;
  auto const d_rcut2 = Utils::sqr(d_rcut);
  auto const d_rcut4 = Utils::sqr(d_rcut2);

  auto const d_a2 = Utils::sqr(alpha_L) / Utils::sqr(box_size);

  auto const d_c = sum_q2 * exp(-d_a2 * d_rcut2);

  d_cc = 4. * Utils::sqr(d_a2) * Utils::sqr(d_rcut2) + 6. * d_a2 * d_rcut2 + 3.;

  d_dc = 8. * Utils::int_pow<3>(d_a2) * Utils::int_pow<3>(d_rcut2) +
         20. * Utils::sqr(d_a2) * d_rcut4 + 30. * d_a2 * d_rcut2 + 15.;

  d_con = 1. / sqrt(Utils::int_pow<3>(box_size) * Utils::sqr(d_a2) * d_rcut *
                    Utils::sqr(d_rcut4) * static_cast<double>(n_c_part));

  d_error_f = d_c * d_con *
              sqrt((13. / 6.) * Utils::sqr(d_cc) +
                   (2. / 15.) * Utils::sqr(d_dc) - (13. / 15.) * d_cc * d_dc);

  return d_error_f;
}

/** Using bisection, find the root of a function "func-tuned_accuracy/sqrt(2.)"
 *  known to lie between x1 and x2. The root, returned as rtbis, will be
 *  refined until its accuracy is \f$\pm\f$ @p xacc.
 */
double dp3m_rtbisection(double box_size, double r_cut_iL, int n_c_part,
                        double sum_q2, double x1, double x2, double xacc,
                        double tuned_accuracy) {
  constexpr int JJ_RTBIS_MAX = 40;

  auto const constant = tuned_accuracy / std::numbers::sqrt2;

  auto const f1 =
      dp3m_real_space_error(box_size, r_cut_iL, n_c_part, sum_q2, x1) -
      constant;
  auto const f2 =
      dp3m_real_space_error(box_size, r_cut_iL, n_c_part, sum_q2, x2) -
      constant;
  if (f1 * f2 >= 0.0) {
    throw std::runtime_error(
        "Root must be bracketed for bisection in dp3m_rtbisection");
  }
  // Orient the search dx, and set rtb to x1 or x2 ...
  double dx;
  double rtb = f1 < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (int j = 1; j <= JJ_RTBIS_MAX; j++) {
    auto const xmid = rtb + (dx *= 0.5);
    auto const fmid =
        dp3m_real_space_error(box_size, r_cut_iL, n_c_part, sum_q2, xmid) -
        constant;
    if (fmid <= 0.0)
      rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0)
      return rtb;
  }
  throw std::runtime_error("Too many bisections in dp3m_rtbisection");
}

void DipolarP3M::sanity_checks_boxl() const {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  for (auto i = 0u; i < 3u; i++) {
    /* check k-space cutoff */
    if (dp3m_params.cao_cut[i] >= box_geo.length_half()[i]) {
      std::stringstream msg;
      msg << "dipolar P3M_init: k-space cutoff " << dp3m_params.cao_cut[i]
          << " is larger than half of box dimension " << box_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
    if (dp3m_params.cao_cut[i] >= local_geo.length()[i]) {
      std::stringstream msg;
      msg << "dipolar P3M_init: k-space cutoff " << dp3m_params.cao_cut[i]
          << " is larger than local box dimension " << local_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
  }

  if ((box_geo.length()[0] != box_geo.length()[1]) or
      (box_geo.length()[1] != box_geo.length()[2])) {
    throw std::runtime_error("DipolarP3M: requires a cubic box");
  }
}

void DipolarP3M::sanity_checks_periodicity() const {
  auto const &box_geo = *get_system().box_geo;
  if (!box_geo.periodic(0) or !box_geo.periodic(1) or !box_geo.periodic(2)) {
    throw std::runtime_error(
        "DipolarP3M: requires periodicity (True, True, True)");
  }
}

void DipolarP3M::sanity_checks_cell_structure() const {
  auto const &local_geo = *get_system().local_geo;
  if (local_geo.cell_structure_type() != CellStructureType::REGULAR and
      local_geo.cell_structure_type() != CellStructureType::HYBRID) {
    throw std::runtime_error(
        "DipolarP3M: requires the regular or hybrid decomposition cell system");
  }
  if (::communicator.size > 1 and
      local_geo.cell_structure_type() == CellStructureType::HYBRID) {
    throw std::runtime_error(
        "DipolarP3M: does not work with the hybrid decomposition cell system, "
        "if using more than one MPI node");
  }
}

void DipolarP3M::sanity_checks_node_grid() const {
  auto const &node_grid = ::communicator.node_grid;
  if (node_grid[0] < node_grid[1] or node_grid[1] < node_grid[2]) {
    throw std::runtime_error(
        "DipolarP3M: node grid must be sorted, largest first");
  }
}

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::scaleby_box_l() {
  auto const &box_geo = *get_system().box_geo;
  dp3m.params.r_cut = dp3m.params.r_cut_iL * box_geo.length()[0];
  dp3m.params.alpha = dp3m.params.alpha_L * box_geo.length_inv()[0];
  dp3m.params.recalc_a_ai_cao_cut(box_geo.length());
  dp3m.local_mesh.recalc_ld_pos(dp3m.params);
  sanity_checks_boxl();
  calc_influence_function_force();
  calc_influence_function_energy();
  dp3m.energy_correction = 0.;
}

template <typename FloatType, Arch Architecture>
void DipolarP3MImpl<FloatType, Architecture>::calc_energy_correction() {
  auto const &box_geo = *get_system().box_geo;
  auto const Ukp3m = calc_average_self_energy_k_space() * box_geo.volume();
  auto const Ewald_volume = Utils::int_pow<3>(dp3m.params.alpha_L);
  auto const Eself = -2. * Ewald_volume * std::numbers::inv_sqrtpi / 3.;
  dp3m.energy_correction =
      -dp3m.sum_mu2 * (Ukp3m + Eself + 2. * std::numbers::pi / 3.);
}

#ifdef NPT
void npt_add_virial_magnetic_contribution(double energy) {
  npt_add_virial_contribution(energy);
}
#endif // NPT

template struct DipolarP3MImpl<float, Arch::CPU>;
template struct DipolarP3MImpl<double, Arch::CPU>;

#endif // DP3M
