/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
 *
 *  The corresponding header file is @ref p3m.hpp.
 */

#include "electrostatics/p3m.hpp"

#ifdef P3M

#include "electrostatics/coulomb.hpp"
#include "electrostatics/elc.hpp"
#include "electrostatics/p3m_gpu.hpp"
#include "electrostatics/p3m_gpu_error.hpp"

#include "p3m/TuningAlgorithm.hpp"
#include "p3m/TuningLogger.hpp"
#include "p3m/fft.hpp"
#include "p3m/influence_function.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "actor/visitors.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "tuning.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/optional.hpp>
#include <boost/range/numeric.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstddef>
#include <functional>
#include <sstream>
#include <stdexcept>

void CoulombP3M::count_charged_particles() {
  auto local_n = 0;
  auto local_q2 = 0.0;
  auto local_q = 0.0;

  for (auto const &p : cell_structure.local_particles()) {
    if (p.q() != 0.0) {
      local_n++;
      local_q2 += Utils::sqr(p.q());
      local_q += p.q();
    }
  }

  boost::mpi::all_reduce(comm_cart, local_n, p3m.sum_qpart, std::plus<>());
  boost::mpi::all_reduce(comm_cart, local_q2, p3m.sum_q2, std::plus<>());
  boost::mpi::all_reduce(comm_cart, local_q, p3m.square_sum_q, std::plus<>());
  p3m.square_sum_q = Utils::sqr(p3m.square_sum_q);
}

/** Calculate the optimal influence function of @cite hockney88a.
 *  (optimised for force calculations)
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft.plan[3].mesh and fft.plan[3].start).
 *
 *  See also: @cite hockney88a eq. 8-22 (p. 275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  @cite deserno98a @cite deserno98b.
 */
void CoulombP3M::calc_influence_function_force() {
  auto const start = Utils::Vector3i{p3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{p3m.fft.plan[3].new_mesh};

  p3m.g_force = grid_influence_function<1>(p3m.params, start, start + size,
                                           box_geo.length());
}

/** Calculate the influence function optimized for the energy and the
 *  self energy correction.
 */
void CoulombP3M::calc_influence_function_energy() {
  auto const start = Utils::Vector3i{p3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{p3m.fft.plan[3].new_mesh};

  p3m.g_energy = grid_influence_function<0>(p3m.params, start, start + size,
                                            box_geo.length());
}

/** Aliasing sum used by @ref p3m_k_space_error. */
static void p3m_tune_aliasing_sums(int nx, int ny, int nz,
                                   Utils::Vector3i const &mesh,
                                   Utils::Vector3d const &mesh_i, int cao,
                                   double alpha_L_i, double *alias1,
                                   double *alias2) {
  using Utils::sinc;

  auto const factor1 = Utils::sqr(Utils::pi() * alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    auto const nmx = nx + mx * mesh[0];
    auto const fnmx = mesh_i[0] * nmx;
    for (int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      auto const nmy = ny + my * mesh[1];
      auto const fnmy = mesh_i[1] * nmy;
      for (int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        auto const nmz = nz + mz * mesh[2];
        auto const fnmz = mesh_i[2] * nmz;

        auto const nm2 = Utils::sqr(nmx) + Utils::sqr(nmy) + Utils::sqr(nmz);
        auto const ex = exp(-factor1 * nm2);
        auto const ex2 = Utils::sqr(ex);

        auto const U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz), 2.0 * cao);

        *alias1 += ex2 / nm2;
        *alias2 += U2 * ex * (nx * nmx + ny * nmy + nz * nmz) / nm2;
      }
    }
  }
}

/** Calculate the real space contribution to the rms error in the force (as
 *  described by Kolafa and Perram).
 *  \param pref       Prefactor of Coulomb interaction.
 *  \param r_cut_iL   rescaled real space cutoff for p3m method.
 *  \param n_c_part   number of charged particles in the system.
 *  \param sum_q2     sum of square of charges in the system
 *  \param alpha_L    rescaled Ewald splitting parameter.
 *  \return real space error
 */
static double p3m_real_space_error(double pref, double r_cut_iL, int n_c_part,
                                   double sum_q2, double alpha_L) {
  return (2. * pref * sum_q2 * exp(-Utils::sqr(r_cut_iL * alpha_L))) /
         sqrt(static_cast<double>(n_c_part) * r_cut_iL * box_geo.length()[0] *
              box_geo.volume());
}

/** Calculate the analytic expression of the error estimate for the
 *  P3M method in @cite hockney88a (eq. 8-23 p. 275) in
 *  order to obtain the rms error in the force for a system of N
 *  randomly distributed particles in a cubic box (k-space part).
 *  \param pref     Prefactor of Coulomb interaction.
 *  \param mesh     number of mesh points in one direction.
 *  \param cao      charge assignment order.
 *  \param n_c_part number of charged particles in the system.
 *  \param sum_q2   sum of square of charges in the system
 *  \param alpha_L  rescaled Ewald splitting parameter.
 *  \return reciprocal (k) space error
 */
static double p3m_k_space_error(double pref, Utils::Vector3i const &mesh,
                                int cao, int n_c_part, double sum_q2,
                                double alpha_L) {
  auto const mesh_i =
      Utils::hadamard_division(Utils::Vector3d::broadcast(1.), mesh);
  auto const alpha_L_i = 1. / alpha_L;
  auto he_q = 0.;

  for (int nx = -mesh[0] / 2; nx < mesh[0] / 2; nx++) {
    auto const ctan_x = p3m_analytic_cotangent_sum(nx, mesh_i[0], cao);
    for (int ny = -mesh[1] / 2; ny < mesh[1] / 2; ny++) {
      auto const ctan_y =
          ctan_x * p3m_analytic_cotangent_sum(ny, mesh_i[1], cao);
      for (int nz = -mesh[2] / 2; nz < mesh[2] / 2; nz++) {
        if ((nx != 0) || (ny != 0) || (nz != 0)) {
          auto const n2 = Utils::sqr(nx) + Utils::sqr(ny) + Utils::sqr(nz);
          auto const cs =
              p3m_analytic_cotangent_sum(nz, mesh_i[2], cao) * ctan_y;
          double alias1, alias2;
          p3m_tune_aliasing_sums(nx, ny, nz, mesh, mesh_i, cao, alpha_L_i,
                                 &alias1, &alias2);

          auto const d = alias1 - Utils::sqr(alias2 / cs) / n2;
          /* at high precision, d can become negative due to extinction;
             also, don't take values that have no significant digits left*/
          if (d > 0 && (fabs(d / alias1) > ROUND_ERROR_PREC))
            he_q += d;
        }
      }
    }
  }
  return 2. * pref * sum_q2 * sqrt(he_q / static_cast<double>(n_c_part)) /
         (box_geo.length()[1] * box_geo.length()[2]);
}

#ifdef CUDA
static double p3mgpu_k_space_error(double prefactor,
                                   Utils::Vector3i const &mesh, int cao,
                                   int npart, double sum_q2, double alpha_L) {
  auto ks_err = 0.;
  if (this_node == 0) {
    ks_err = p3m_k_space_error_gpu(prefactor, mesh.data(), cao, npart, sum_q2,
                                   alpha_L, box_geo.length().data());
  }
  boost::mpi::broadcast(comm_cart, ks_err, 0);
  return ks_err;
}
#endif

void CoulombP3M::init() {
  assert(p3m.params.mesh >= Utils::Vector3i::broadcast(1));
  assert(p3m.params.cao >= 1 and p3m.params.cao <= 7);
  assert(p3m.params.alpha > 0.);

  p3m.params.cao3 = Utils::int_pow<3>(p3m.params.cao);
  p3m.params.recalc_a_ai_cao_cut(box_geo.length());

  sanity_checks();

  double elc_layer = 0.;
  if (auto elc_actor = get_actor_by_type<ElectrostaticLayerCorrection>(
          electrostatics_actor)) {
    elc_layer = elc_actor->elc.space_layer;
  }

  p3m.local_mesh.calc_local_ca_mesh(p3m.params, local_geo, skin, elc_layer);
  p3m.sm.resize(comm_cart, p3m.local_mesh);

  int ca_mesh_size =
      fft_init(p3m.local_mesh.dim, p3m.local_mesh.margin, p3m.params.mesh,
               p3m.params.mesh_off, p3m.ks_pnum, p3m.fft, node_grid, comm_cart);
  p3m.rs_mesh.resize(ca_mesh_size);

  for (auto &e : p3m.E_mesh) {
    e.resize(ca_mesh_size);
  }

  p3m.calc_differential_operator();

  /* fix box length dependent constants */
  scaleby_box_l();

  count_charged_particles();
}

CoulombP3M::CoulombP3M(P3MParameters &&parameters, double prefactor,
                       int tune_timings, bool tune_verbose,
                       bool check_complex_residuals)
    : p3m{std::move(parameters)}, tune_timings{tune_timings},
      tune_verbose{tune_verbose}, check_complex_residuals{
                                      check_complex_residuals} {

  if (tune_timings <= 0) {
    throw std::domain_error("Parameter 'timings' must be > 0");
  }
  m_is_tuned = !p3m.params.tuning;
  p3m.params.tuning = false;
  set_prefactor(prefactor);
}

namespace {
template <int cao> struct AssignCharge {
  void operator()(p3m_data_struct &p3m, double q,
                  Utils::Vector3d const &real_pos,
                  p3m_interpolation_cache &inter_weights) {
    auto const w = p3m_calculate_interpolation_weights<cao>(
        real_pos, p3m.params.ai, p3m.local_mesh);

    inter_weights.store(w);

    p3m_interpolate(p3m.local_mesh, w, [q, &p3m](int ind, double w) {
      p3m.rs_mesh[ind] += w * q;
    });
  }

  void operator()(p3m_data_struct &p3m, double q,
                  Utils::Vector3d const &real_pos) {
    p3m_interpolate(
        p3m.local_mesh,
        p3m_calculate_interpolation_weights<cao>(real_pos, p3m.params.ai,
                                                 p3m.local_mesh),
        [q, &p3m](int ind, double w) { p3m.rs_mesh[ind] += w * q; });
  }

  void operator()(p3m_data_struct &p3m, ParticleRange const &particles) {
    for (auto &p : particles) {
      if (p.q() != 0.0) {
        this->operator()(p3m, p.q(), p.pos(), p3m.inter_weights);
      }
    }
  }
};
} // namespace

void CoulombP3M::charge_assign(ParticleRange const &particles) {
  p3m.inter_weights.reset(p3m.params.cao);

  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  Utils::integral_parameter<int, AssignCharge, 1, 7>(p3m.params.cao, p3m,
                                                     particles);
}

void CoulombP3M::assign_charge(double q, Utils::Vector3d const &real_pos,
                               p3m_interpolation_cache &inter_weights) {
  Utils::integral_parameter<int, AssignCharge, 1, 7>(p3m.params.cao, p3m, q,
                                                     real_pos, inter_weights);
}

void CoulombP3M::assign_charge(double q, Utils::Vector3d const &real_pos) {
  Utils::integral_parameter<int, AssignCharge, 1, 7>(p3m.params.cao, p3m, q,
                                                     real_pos);
}

namespace {
template <int cao> struct AssignForces {
  void operator()(p3m_data_struct &p3m, double force_prefac,
                  ParticleRange const &particles) const {
    using Utils::make_const_span;
    using Utils::Span;
    using Utils::Vector;

    assert(cao == p3m.inter_weights.cao());

    /* charged particle counter */
    auto p_index = std::size_t{0ul};

    for (auto &p : particles) {
      if (p.q() != 0.0) {
        auto const pref = p.q() * force_prefac;
        auto const w = p3m.inter_weights.load<cao>(p_index);

        Utils::Vector3d force{};
        p3m_interpolate(p3m.local_mesh, w, [&force, &p3m](int ind, double w) {
          force += w * Utils::Vector3d{p3m.E_mesh[0][ind], p3m.E_mesh[1][ind],
                                       p3m.E_mesh[2][ind]};
        });

        p.force() -= pref * force;
        ++p_index;
      }
    }
  }
};

auto dipole_moment(Particle const &p, BoxGeometry const &box) {
  return p.q() * unfolded_position(p.pos(), p.image_box(), box.length());
}

auto calc_dipole_moment(boost::mpi::communicator const &comm,
                        ParticleRange const &particles,
                        BoxGeometry const &box) {
  auto const local_dip = boost::accumulate(
      particles, Utils::Vector3d{}, [&box](Utils::Vector3d dip, auto const &p) {
        return dip + dipole_moment(p, box);
      });

  return boost::mpi::all_reduce(comm, local_dip, std::plus<>());
}
} // namespace

/** @details Calculate the long range electrostatics part of the pressure
 *  tensor. This is part \f$\Pi_{\textrm{dir}, \alpha, \beta}\f$ eq. (2.6)
 *  in @cite essmann95a. The part \f$\Pi_{\textrm{corr}, \alpha, \beta}\f$
 *  eq. (2.8) is not present here since M is the empty set in our simulations.
 */
Utils::Vector9d CoulombP3M::p3m_calc_kspace_pressure_tensor() {
  using namespace detail::FFT_indexing;

  Utils::Vector9d node_k_space_pressure_tensor{};

  if (p3m.sum_q2 > 0.) {
    p3m.sm.gather_grid(p3m.rs_mesh.data(), comm_cart, p3m.local_mesh.dim);
    fft_perform_forw(p3m.rs_mesh.data(), p3m.fft, comm_cart);

    auto diagonal = 0.;
    int ind = 0;
    int j[3];
    auto const half_alpha_inv_sq = Utils::sqr(1. / 2. / p3m.params.alpha);
    for (j[0] = 0; j[0] < p3m.fft.plan[3].new_mesh[RX]; j[0]++) {
      for (j[1] = 0; j[1] < p3m.fft.plan[3].new_mesh[RY]; j[1]++) {
        for (j[2] = 0; j[2] < p3m.fft.plan[3].new_mesh[RZ]; j[2]++) {
          auto const kx = 2. * Utils::pi() *
                          p3m.d_op[RX][j[KX] + p3m.fft.plan[3].start[KX]] *
                          box_geo.length_inv()[RX];
          auto const ky = 2. * Utils::pi() *
                          p3m.d_op[RY][j[KY] + p3m.fft.plan[3].start[KY]] *
                          box_geo.length_inv()[RY];
          auto const kz = 2. * Utils::pi() *
                          p3m.d_op[RZ][j[KZ] + p3m.fft.plan[3].start[KZ]] *
                          box_geo.length_inv()[RZ];
          auto const sqk = Utils::sqr(kx) + Utils::sqr(ky) + Utils::sqr(kz);

          if (sqk != 0.) {
            auto const node_k_space_energy =
                p3m.g_energy[ind] * (Utils::sqr(p3m.rs_mesh[2 * ind]) +
                                     Utils::sqr(p3m.rs_mesh[2 * ind + 1]));
            auto const vterm = -2. * (1. / sqk + half_alpha_inv_sq);
            auto const pref = node_k_space_energy * vterm;
            node_k_space_pressure_tensor[0] += pref * kx * kx; /* sigma_xx */
            node_k_space_pressure_tensor[1] += pref * kx * ky; /* sigma_xy */
            node_k_space_pressure_tensor[2] += pref * kx * kz; /* sigma_xz */
            node_k_space_pressure_tensor[3] += pref * ky * kx; /* sigma_yx */
            node_k_space_pressure_tensor[4] += pref * ky * ky; /* sigma_yy */
            node_k_space_pressure_tensor[5] += pref * ky * kz; /* sigma_yz */
            node_k_space_pressure_tensor[6] += pref * kz * kx; /* sigma_zx */
            node_k_space_pressure_tensor[7] += pref * kz * ky; /* sigma_zy */
            node_k_space_pressure_tensor[8] += pref * kz * kz; /* sigma_zz */
            diagonal += node_k_space_energy;
          }
          ind++;
        }
      }
    }
    node_k_space_pressure_tensor[0] += diagonal;
    node_k_space_pressure_tensor[4] += diagonal;
    node_k_space_pressure_tensor[8] += diagonal;
  }

  return node_k_space_pressure_tensor * prefactor / (2. * box_geo.volume());
}

double CoulombP3M::long_range_kernel(bool force_flag, bool energy_flag,
                                     ParticleRange const &particles) {
  /* Gather information for FFT grid inside the nodes domain (inner local mesh)
   * and perform forward 3D FFT (Charge Assignment Mesh). */
  p3m.sm.gather_grid(p3m.rs_mesh.data(), comm_cart, p3m.local_mesh.dim);
  fft_perform_forw(p3m.rs_mesh.data(), p3m.fft, comm_cart);

  // Note: after these calls, the grids are in the order yzx and not xyz
  // anymore!!!
  /* The dipole moment is only needed if we don't have metallic boundaries. */
  auto const box_dipole = (p3m.params.epsilon != P3M_EPSILON_METALLIC)
                              ? boost::make_optional(calc_dipole_moment(
                                    comm_cart, particles, box_geo))
                              : boost::none;
  auto const volume = box_geo.volume();
  auto const pref = 4. * Utils::pi() / volume / (2. * p3m.params.epsilon + 1.);

  /* === k-space force calculation  === */
  if (force_flag) {
    /* sqrt(-1)*k differentiation */
    int j[3];
    int ind = 0;
    for (j[0] = 0; j[0] < p3m.fft.plan[3].new_mesh[0]; j[0]++) {
      for (j[1] = 0; j[1] < p3m.fft.plan[3].new_mesh[1]; j[1]++) {
        for (j[2] = 0; j[2] < p3m.fft.plan[3].new_mesh[2]; j[2]++) {
          auto const rho_hat = std::complex<double>(p3m.rs_mesh[2 * ind + 0],
                                                    p3m.rs_mesh[2 * ind + 1]);
          auto const phi_hat = p3m.g_force[ind] * rho_hat;

          for (int d = 0; d < 3; d++) {
            /* direction in r-space: */
            int d_rs = (d + p3m.ks_pnum) % 3;
            /* directions */
            auto const k = 2. * Utils::pi() *
                           p3m.d_op[d_rs][j[d] + p3m.fft.plan[3].start[d]] *
                           box_geo.length_inv()[d_rs];

            /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
            p3m.E_mesh[d_rs][2 * ind + 0] = -k * phi_hat.imag();
            p3m.E_mesh[d_rs][2 * ind + 1] = +k * phi_hat.real();
          }

          ind++;
        }
      }
    }

    /* Back FFT force component mesh */
    auto const check_complex = !p3m.params.tuning and check_complex_residuals;
    for (int d = 0; d < 3; d++) {
      fft_perform_back(p3m.E_mesh[d].data(), check_complex, p3m.fft, comm_cart);
    }

    /* redistribute force component mesh */
    std::array<double *, 3> E_fields = {
        {p3m.E_mesh[0].data(), p3m.E_mesh[1].data(), p3m.E_mesh[2].data()}};
    p3m.sm.spread_grid(Utils::make_span(E_fields), comm_cart,
                       p3m.local_mesh.dim);

    auto const force_prefac = prefactor / volume;
    Utils::integral_parameter<int, AssignForces, 1, 7>(p3m.params.cao, p3m,
                                                       force_prefac, particles);

    // add dipole forces
    if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
      auto const dm = prefactor * pref * box_dipole.value();
      for (auto &p : particles) {
        p.force() -= p.q() * dm;
      }
    }
  }

  /* === k-space energy calculation  === */
  if (energy_flag) {
    auto node_energy = 0.;
    for (int i = 0; i < p3m.fft.plan[3].new_size; i++) {
      // Use the energy optimized influence function for energy!
      node_energy += p3m.g_energy[i] * (Utils::sqr(p3m.rs_mesh[2 * i]) +
                                        Utils::sqr(p3m.rs_mesh[2 * i + 1]));
    }
    node_energy /= 2. * volume;

    auto energy = 0.;
    boost::mpi::reduce(comm_cart, node_energy, energy, std::plus<>(), 0);
    if (this_node == 0) {
      /* self energy correction */
      energy -= p3m.sum_q2 * p3m.params.alpha * Utils::sqrt_pi_i();
      /* net charge correction */
      energy -= p3m.square_sum_q * Utils::pi() /
                (2. * volume * Utils::sqr(p3m.params.alpha));
      /* dipole correction */
      if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
        energy += pref * box_dipole.value().norm2();
      }
    }
    return prefactor * energy;
  }

  return 0.;
}

class CoulombTuningAlgorithm : public TuningAlgorithm {
  p3m_data_struct &p3m;
  double m_mesh_density_min = -1., m_mesh_density_max = -1.;
  // indicates if mesh should be tuned
  bool m_tune_mesh = false;

public:
  CoulombTuningAlgorithm(p3m_data_struct &input_p3m, double prefactor,
                         int timings)
      : TuningAlgorithm{prefactor, timings}, p3m{input_p3m} {}

  P3MParameters &get_params() override { return p3m.params; }

  void on_solver_change() const override { on_coulomb_change(); }

  void setup_logger(bool verbose) override {
#ifdef CUDA
    auto const on_gpu = has_actor_of_type<CoulombP3MGPU>(electrostatics_actor);
#else
    auto const on_gpu = false;
#endif
    m_logger = std::make_unique<TuningLogger>(
        verbose and this_node == 0, (on_gpu) ? "CoulombP3MGPU" : "CoulombP3M",
        TuningLogger::Mode::Coulomb);
    m_logger->tuning_goals(p3m.params.accuracy, m_prefactor,
                           box_geo.length()[0], p3m.sum_qpart, p3m.sum_q2);
    m_logger->log_tuning_start();
  }

  boost::optional<std::string>
  layer_correction_veto_r_cut(double r_cut) const override {
    if (auto elc_actor = get_actor_by_type<ElectrostaticLayerCorrection>(
            electrostatics_actor)) {
      return elc_actor->veto_r_cut(r_cut);
    }
    return {};
  }

  std::tuple<double, double, double, double>
  calculate_accuracy(Utils::Vector3i const &mesh, int cao,
                     double r_cut_iL) const override {

    double alpha_L, rs_err, ks_err;

    /* calc maximal real space error for setting */
    rs_err = p3m_real_space_error(m_prefactor, r_cut_iL, p3m.sum_qpart,
                                  p3m.sum_q2, 0.);

    if (Utils::sqrt_2() * rs_err > p3m.params.accuracy) {
      /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
      alpha_L =
          sqrt(log(Utils::sqrt_2() * rs_err / p3m.params.accuracy)) / r_cut_iL;
    } else {
      /* even alpha=0 is ok, however, we cannot choose it since it kills the
         k-space error formula.
         Anyways, this very likely NOT the optimal solution */
      alpha_L = 0.1;
    }

    /* calculate real-space and k-space error for this alpha_L */
    rs_err = p3m_real_space_error(m_prefactor, r_cut_iL, p3m.sum_qpart,
                                  p3m.sum_q2, alpha_L);
#ifdef CUDA
    if (has_actor_of_type<CoulombP3MGPU>(electrostatics_actor)) {
      ks_err = p3mgpu_k_space_error(m_prefactor, mesh, cao, p3m.sum_qpart,
                                    p3m.sum_q2, alpha_L);
    } else
#endif
      ks_err = p3m_k_space_error(m_prefactor, mesh, cao, p3m.sum_qpart,
                                 p3m.sum_q2, alpha_L);

    return {Utils::Vector2d{rs_err, ks_err}.norm(), rs_err, ks_err, alpha_L};
  }

  void determine_mesh_limits() override {
    auto const mesh_density =
        static_cast<double>(p3m.params.mesh[0]) * box_geo.length_inv()[0];

    if (p3m.params.mesh == Utils::Vector3i::broadcast(-1)) {
      /* avoid using more than 1 GB of FFT arrays */
      auto const normalized_box_dim = std::cbrt(box_geo.volume());
      auto const max_npart_per_dim = 512.;
      /* simple heuristic to limit the tried meshes if the accuracy cannot
         be obtained with smaller meshes, but normally not all these
         meshes have to be tested */
      auto const min_npart_per_dim = std::min(
          max_npart_per_dim, std::cbrt(static_cast<double>(p3m.sum_qpart)));
      m_mesh_density_min = min_npart_per_dim / normalized_box_dim;
      m_mesh_density_max = max_npart_per_dim / normalized_box_dim;
      m_tune_mesh = true;
    } else {
      m_mesh_density_min = m_mesh_density_max = mesh_density;
      assert(p3m.params.mesh[0] >= 1);
      if (p3m.params.mesh[1] == -1 and p3m.params.mesh[2] == -1) {
        // determine the two missing values by rescaling by the box length
        for (int i : {1, 2}) {
          p3m.params.mesh[i] =
              static_cast<int>(std::round(mesh_density * box_geo.length()[i]));
          // make the mesh even in all directions
          p3m.params.mesh[i] += p3m.params.mesh[i] % 2;
        }
      }
      m_logger->report_fixed_mesh(p3m.params.mesh);
    }
  }

  TuningAlgorithm::Parameters get_time() override {
    auto tuned_params = TuningAlgorithm::Parameters{};
    auto time_best = time_sentinel;
    auto mesh_density = m_mesh_density_min;
    while (mesh_density <= m_mesh_density_max) {
      auto trial_params = TuningAlgorithm::Parameters{};
      if (m_tune_mesh) {
        for (int i : {0, 1, 2}) {
          trial_params.mesh[i] =
              static_cast<int>(std::round(box_geo.length()[i] * mesh_density));
          // make the mesh even in all directions
          trial_params.mesh[i] += trial_params.mesh[i] % 2;
        }
      } else {
        trial_params.mesh = p3m.params.mesh;
      }
      trial_params.cao = cao_best;

      auto const trial_time =
          get_m_time(trial_params.mesh, trial_params.cao, trial_params.r_cut_iL,
                     trial_params.alpha_L, trial_params.accuracy);

      if (trial_time >= 0.) {
        /* the optimum r_cut for this mesh is the upper limit for higher meshes,
           everything else is slower */
        if (has_actor_of_type<CoulombP3M>(electrostatics_actor)) {
          m_r_cut_iL_max = trial_params.r_cut_iL;
        }

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
      mesh_density += 0.1;
    }
    return tuned_params;
  }
};

void CoulombP3M::tune() {
  if (p3m.params.alpha_L == 0. and p3m.params.alpha != 0.) {
    p3m.params.alpha_L = p3m.params.alpha * box_geo.length()[0];
  }
  if (p3m.params.r_cut_iL == 0. and p3m.params.r_cut != 0.) {
    p3m.params.r_cut_iL = p3m.params.r_cut * box_geo.length_inv()[0];
  }
  if (not is_tuned()) {
    count_charged_particles();
    if (p3m.sum_qpart == 0) {
      throw std::runtime_error(
          "CoulombP3M: no charged particles in the system");
    }
    try {
      CoulombTuningAlgorithm parameters(p3m, prefactor, tune_timings);
      parameters.setup_logger(tune_verbose);
      // parameter ranges
      parameters.determine_mesh_limits();
      parameters.determine_r_cut_limits();
      parameters.determine_cao_limits(7);
      // run tuning algorithm
      parameters.tune();
      m_is_tuned = true;
      on_coulomb_change();
    } catch (...) {
      p3m.params.tuning = false;
      throw;
    }
  }
  init();
}

void CoulombP3M::sanity_checks_boxl() const {
  for (unsigned int i = 0; i < 3; i++) {
    /* check k-space cutoff */
    if (p3m.params.cao_cut[i] >= box_geo.length_half()[i]) {
      std::stringstream msg;
      msg << "P3M_init: k-space cutoff " << p3m.params.cao_cut[i]
          << " is larger than half of box dimension " << box_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
    if (p3m.params.cao_cut[i] >= local_geo.length()[i]) {
      std::stringstream msg;
      msg << "P3M_init: k-space cutoff " << p3m.params.cao_cut[i]
          << " is larger than local box dimension " << local_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
  }

  if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
    if ((box_geo.length()[0] != box_geo.length()[1]) or
        (box_geo.length()[1] != box_geo.length()[2]) or
        (p3m.params.mesh[0] != p3m.params.mesh[1]) or
        (p3m.params.mesh[1] != p3m.params.mesh[2])) {
      throw std::runtime_error(
          "CoulombP3M: non-metallic epsilon requires cubic box");
    }
  }
}

void CoulombP3M::sanity_checks_periodicity() const {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error(
        "CoulombP3M: requires periodicity (True, True, True)");
  }
}

void CoulombP3M::sanity_checks_cell_structure() const {
  if (local_geo.cell_structure_type() !=
          CellStructureType::CELL_STRUCTURE_REGULAR &&
      local_geo.cell_structure_type() !=
          CellStructureType::CELL_STRUCTURE_HYBRID) {
    throw std::runtime_error(
        "CoulombP3M: requires the regular or hybrid decomposition cell system");
  }
  if (n_nodes > 1 && local_geo.cell_structure_type() ==
                         CellStructureType::CELL_STRUCTURE_HYBRID) {
    throw std::runtime_error(
        "CoulombP3M: does not work with the hybrid decomposition cell system, "
        "if using more than one MPI node");
  }
}

void CoulombP3M::sanity_checks_node_grid() const {
  if (node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    throw std::runtime_error(
        "CoulombP3M: node grid must be sorted, largest first");
  }
}

void CoulombP3M::scaleby_box_l() {
  p3m.params.r_cut = p3m.params.r_cut_iL * box_geo.length()[0];
  p3m.params.alpha = p3m.params.alpha_L * box_geo.length_inv()[0];
  p3m.params.recalc_a_ai_cao_cut(box_geo.length());
  p3m.local_mesh.recalc_ld_pos(p3m.params);
  sanity_checks_boxl();
  calc_influence_function_force();
  calc_influence_function_energy();
}

#endif // P3M
