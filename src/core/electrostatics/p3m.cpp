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
 *
 *  The corresponding header file is @ref p3m.hpp.
 */

#include "config/config.hpp"

#ifdef P3M

#include "electrostatics/p3m.hpp"
#include "electrostatics/p3m.impl.hpp"

#include "electrostatics/coulomb.hpp"
#include "electrostatics/elc.hpp"
#ifdef CUDA
#include "electrostatics/p3m_gpu_cuda.cuh"
#include "electrostatics/p3m_gpu_error.hpp"
#endif // CUDA

#include "fft/fft.hpp"
#include "p3m/TuningAlgorithm.hpp"
#include "p3m/TuningLogger.hpp"
#include "p3m/for_each_3d.hpp"
#include "p3m/influence_function.hpp"
#include "p3m/math.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticlePropertyIterator.hpp"
#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "actor/visitors.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrators/Propagation.hpp"
#include "npt.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <numbers>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#ifdef FFTW3_H
#error "The FFTW3 library shouldn't be visible in this translation unit"
#endif

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::count_charged_particles() {
  auto local_n = 0;
  auto local_q2 = 0.0;
  auto local_q = 0.0;

  for (auto const &p : get_system().cell_structure->local_particles()) {
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
 *  Each node calculates only the values for its domain in k-space.
 *
 *  See also: @cite hockney88a eq. 8-22 (p. 275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  @cite deserno98a @cite deserno98b.
 */
template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::calc_influence_function_force() {
  auto const [KX, KY, KZ] = p3m.fft->get_permutations();
  p3m.g_force = grid_influence_function<FloatType, 1>(
      p3m.params, p3m.mesh.start, p3m.mesh.stop, KX, KY, KZ,
      get_system().box_geo->length_inv());
}

/** Calculate the influence function optimized for the energy and the
 *  self energy correction.
 */
template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::calc_influence_function_energy() {
  auto const [KX, KY, KZ] = p3m.fft->get_permutations();
  p3m.g_energy = grid_influence_function<FloatType, 0>(
      p3m.params, p3m.mesh.start, p3m.mesh.stop, KX, KY, KZ,
      get_system().box_geo->length_inv());
}

/** Aliasing sum used by @ref p3m_k_space_error. */
static auto p3m_tune_aliasing_sums(Utils::Vector3i const &shift,
                                   Utils::Vector3i const &mesh,
                                   Utils::Vector3d const &mesh_i, int cao,
                                   double alpha_L_i) {

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
        auto const ex = exp(-factor1 * norm_sq);
        auto const energy = std::pow(Utils::product(fnm), 2 * cao);
        alias1 += Utils::sqr(ex) / norm_sq;
        alias2 += energy * ex * (shift * nm) / norm_sq;
      },
      [&](unsigned dim, int n) {
        nm[dim] = shift[dim] + n * mesh[dim];
        fnm[dim] = math::sinc(nm[dim] * mesh_i[dim]);
      });

  return std::make_pair(alias1, alias2);
}

/** Calculate the real space contribution to the rms error in the force (as
 *  described by Kolafa and Perram).
 *  \param pref       Prefactor of Coulomb interaction.
 *  \param r_cut_iL   rescaled real space cutoff for p3m method.
 *  \param n_c_part   number of charged particles in the system.
 *  \param sum_q2     sum of square of charges in the system
 *  \param alpha_L    rescaled Ewald splitting parameter.
 *  \param box_l      box dimensions.
 *  \return real space error
 */
static double p3m_real_space_error(double pref, double r_cut_iL, int n_c_part,
                                   double sum_q2, double alpha_L,
                                   Utils::Vector3d const &box_l) {
  auto const volume = Utils::product(box_l);
  return (2. * pref * sum_q2 * exp(-Utils::sqr(r_cut_iL * alpha_L))) /
         sqrt(static_cast<double>(n_c_part) * r_cut_iL * box_l[0] * volume);
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
 *  \param box_l    box dimensions.
 *  \return reciprocal (k) space error
 */
static double p3m_k_space_error(double pref, Utils::Vector3i const &mesh,
                                int cao, int n_c_part, double sum_q2,
                                double alpha_L, Utils::Vector3d const &box_l) {

  auto const cotangent_sum = math::get_analytic_cotangent_sum_kernel(cao);
  auto const mesh_i = 1. / Utils::Vector3d(mesh);
  auto const alpha_L_i = 1. / alpha_L;
  auto const mesh_stop = mesh / 2;
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
              p3m_tune_aliasing_sums(indices, mesh, mesh_i, cao, alpha_L_i);
          auto const d = alias1 - Utils::sqr(alias2 / cs) / n2;
          /* at high precision, d can become negative due to extinction;
             also, don't take values that have no significant digits left*/
          if (d > 0. and std::fabs(d / alias1) > ROUND_ERROR_PREC) {
            he_q += d;
          }
        }
      },
      [&values, &mesh_i, cotangent_sum](unsigned dim, int n) {
        values[dim] = cotangent_sum(n, mesh_i[dim]);
      });

  return 2. * pref * sum_q2 * sqrt(he_q / static_cast<double>(n_c_part)) /
         (box_l[1] * box_l[2]);
}

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::init_cpu_kernels() {
  assert(p3m.params.mesh >= Utils::Vector3i::broadcast(1));
  assert(p3m.params.cao >= 1 and p3m.params.cao <= 7);
  assert(p3m.params.alpha > 0.);

  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  auto const skin = system.cell_structure->get_verlet_skin();

  p3m.params.cao3 = Utils::int_pow<3>(p3m.params.cao);
  p3m.params.recalc_a_ai_cao_cut(box_geo.length());

  sanity_checks();

  auto const &solver = system.coulomb.impl->solver;
  double elc_layer = 0.;
  if (auto actor = get_actor_by_type<ElectrostaticLayerCorrection>(solver)) {
    elc_layer = actor->elc.space_layer;
  }

  assert(p3m.fft);
  p3m.local_mesh.calc_local_ca_mesh(p3m.params, local_geo, skin, elc_layer);
  p3m.fft_buffers->init_halo();
  p3m.fft->init(p3m.params);
  p3m.mesh.ks_pnum = p3m.fft->get_ks_pnum();
  p3m.fft_buffers->init_meshes(p3m.fft->get_ca_mesh_size());
  p3m.update_mesh_views();
  p3m.calc_differential_operator();

  /* fix box length dependent constants */
  scaleby_box_l();

  count_charged_particles();
}

namespace {
template <int cao> struct AssignCharge {
  void operator()(auto &p3m, double q, Utils::Vector3d const &real_pos,
                  InterpolationWeights<cao> const &w) {
    using value_type =
        typename std::remove_reference_t<decltype(p3m)>::value_type;
    p3m_interpolate(p3m.local_mesh, w, [q, &p3m](int ind, double w) {
      p3m.mesh.rs_scalar[ind] += value_type(w * q);
    });
  }

  void operator()(auto &p3m, double q, Utils::Vector3d const &real_pos,
                  p3m_interpolation_cache &inter_weights) {
    auto const w = p3m_calculate_interpolation_weights<cao>(
        real_pos, p3m.params.ai, p3m.local_mesh);
    inter_weights.store(w);
    this->operator()(p3m, q, real_pos, w);
  }

  void operator()(auto &p3m, double q, Utils::Vector3d const &real_pos) {
    auto const w = p3m_calculate_interpolation_weights<cao>(
        real_pos, p3m.params.ai, p3m.local_mesh);
    this->operator()(p3m, q, real_pos, w);
  }

  template <typename combined_ranges>
  void operator()(auto &p3m, combined_ranges const &p_q_pos_range) {
    for (auto zipped : p_q_pos_range) {
      auto const p_q = boost::get<0>(zipped);
      auto const &p_pos = boost::get<1>(zipped);
      if (p_q != 0.0) {
        this->operator()(p3m, p_q, p_pos, p3m.inter_weights);
      }
    }
  }
};
} // namespace

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::charge_assign(
    ParticleRange const &particles) {
  prepare_fft_mesh(true);

  auto p_q_range = ParticlePropertyRange::charge_range(particles);
  auto p_pos_range = ParticlePropertyRange::pos_range(particles);

  Utils::integral_parameter<int, AssignCharge, 1, 7>(
      p3m.params.cao, p3m, boost::combine(p_q_range, p_pos_range));
}

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::assign_charge(
    double q, Utils::Vector3d const &real_pos, bool skip_cache) {
  if (skip_cache) {
    Utils::integral_parameter<int, AssignCharge, 1, 7>(p3m.params.cao, p3m, q,
                                                       real_pos);
  } else {
    Utils::integral_parameter<int, AssignCharge, 1, 7>(
        p3m.params.cao, p3m, q, real_pos, p3m.inter_weights);
  }
}

template <int cao> struct AssignForces {
  template <typename combined_ranges>
  void operator()(auto &p3m, double force_prefac,
                  combined_ranges const &p_q_force_range) const {

    assert(cao == p3m.inter_weights.cao());

    /* charged particle counter */
    auto p_index = std::size_t{0ul};

    for (auto zipped : p_q_force_range) {
      auto p_q = boost::get<0>(zipped);
      auto &p_force = boost::get<1>(zipped);
      if (p_q != 0.0) {
        auto const pref = p_q * force_prefac;
        auto const w = p3m.inter_weights.template load<cao>(p_index);

        Utils::Vector3d force{};
        p3m_interpolate(p3m.local_mesh, w, [&force, &p3m](int ind, double w) {
          force[0u] += w * double(p3m.mesh.rs_fields[0u][ind]);
          force[1u] += w * double(p3m.mesh.rs_fields[1u][ind]);
          force[2u] += w * double(p3m.mesh.rs_fields[2u][ind]);
        });

        p_force -= pref * force;
        ++p_index;
      }
    }
  }
};

template <typename combined_ranges>
static auto calc_dipole_moment(boost::mpi::communicator const &comm,
                               combined_ranges const &p_q_unfolded_pos_range) {
  auto const local_dip =
      boost::accumulate(p_q_unfolded_pos_range, Utils::Vector3d{},
                        [](Utils::Vector3d const &dip, auto const &q_pos) {
                          auto const p_q = boost::get<0>(q_pos);
                          auto const &p_unfolded_pos = boost::get<1>(q_pos);
                          return dip + p_q * p_unfolded_pos;
                        });
  return boost::mpi::all_reduce(comm, local_dip, std::plus<>());
}

/** @details Calculate the long range electrostatics part of the pressure
 *  tensor. This is part \f$\Pi_{\textrm{dir}, \alpha, \beta}\f$ eq. (2.6)
 *  in @cite essmann95a. The part \f$\Pi_{\textrm{corr}, \alpha, \beta}\f$
 *  eq. (2.8) is not present here since M is the empty set in our simulations.
 */
template <typename FloatType, Arch Architecture>
Utils::Vector9d CoulombP3MImpl<FloatType, Architecture>::long_range_pressure(
    ParticleRange const &particles) {
  auto const &box_geo = *get_system().box_geo;
  Utils::Vector9d node_k_space_pressure_tensor{};

  if (p3m.sum_q2 > 0.) {
    charge_assign(particles);
    p3m.fft_buffers->perform_scalar_halo_gather();
    p3m.fft->forward_fft(p3m.fft_buffers->get_scalar_mesh());
    p3m.update_mesh_views();

    auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
    auto const &mesh_stop = p3m.mesh.size;
    auto const &offset = p3m.mesh.start;
    auto const half_alpha_inv_sq = Utils::sqr(1. / 2. / p3m.params.alpha);
    auto const wavevector = (2. * std::numbers::pi) * box_geo.length_inv();
    auto const [KX, KY, KZ] = p3m.fft->get_permutations();
    auto indices = Utils::Vector3i{};
    auto index = std::size_t(0u);
    auto diagonal = 0.;

    for_each_3d(mesh_start, mesh_stop, indices, [&]() {
      auto const shift = indices + offset;
      auto const kx = p3m.d_op[0u][shift[KX]] * wavevector[0u];
      auto const ky = p3m.d_op[1u][shift[KY]] * wavevector[1u];
      auto const kz = p3m.d_op[2u][shift[KZ]] * wavevector[2u];
      auto const norm_sq = Utils::sqr(kx) + Utils::sqr(ky) + Utils::sqr(kz);

      if (norm_sq != 0.) {
        auto const node_k_space_energy =
            double(p3m.g_energy[index] *
                   (Utils::sqr(p3m.mesh.rs_scalar[2u * index + 0u]) +
                    Utils::sqr(p3m.mesh.rs_scalar[2u * index + 1u])));
        auto const vterm = -2. * (1. / norm_sq + half_alpha_inv_sq);
        auto const pref = node_k_space_energy * vterm;
        node_k_space_pressure_tensor[0u] += pref * kx * kx; /* sigma_xx */
        node_k_space_pressure_tensor[1u] += pref * kx * ky; /* sigma_xy */
        node_k_space_pressure_tensor[2u] += pref * kx * kz; /* sigma_xz */
        node_k_space_pressure_tensor[4u] += pref * ky * ky; /* sigma_yy */
        node_k_space_pressure_tensor[5u] += pref * ky * kz; /* sigma_yz */
        node_k_space_pressure_tensor[8u] += pref * kz * kz; /* sigma_zz */
        diagonal += node_k_space_energy;
      }
      ++index;
    });

    node_k_space_pressure_tensor[0u] += diagonal;
    node_k_space_pressure_tensor[4u] += diagonal;
    node_k_space_pressure_tensor[8u] += diagonal;
    node_k_space_pressure_tensor[3u] = node_k_space_pressure_tensor[1u];
    node_k_space_pressure_tensor[6u] = node_k_space_pressure_tensor[2u];
    node_k_space_pressure_tensor[7u] = node_k_space_pressure_tensor[5u];
  }

  return node_k_space_pressure_tensor * prefactor / (2. * box_geo.volume());
}

template <typename FloatType, Arch Architecture>
double CoulombP3MImpl<FloatType, Architecture>::long_range_kernel(
    bool force_flag, bool energy_flag, ParticleRange const &particles) {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
#ifdef NPT
  auto const npt_flag =
      force_flag and (system.propagation->integ_switch == INTEG_METHOD_NPT_ISO);
#else
  auto constexpr npt_flag = false;
#endif

  if (p3m.sum_q2 > 0.) {
    if (not has_actor_of_type<ElectrostaticLayerCorrection>(
            system.coulomb.impl->solver)) {
      charge_assign(particles);
    }
    p3m.fft_buffers->perform_scalar_halo_gather();
    p3m.fft->forward_fft(p3m.fft_buffers->get_scalar_mesh());
    p3m.update_mesh_views();
  }

  auto p_q_range = ParticlePropertyRange::charge_range(particles);
  auto p_force_range = ParticlePropertyRange::force_range(particles);
  auto p_unfolded_pos_range =
      ParticlePropertyRange::unfolded_pos_range(particles, box_geo);

  // Note: after these calls, the grids are in the order yzx and not xyz
  // anymore!!!
  /* The dipole moment is only needed if we don't have metallic boundaries. */
  auto const box_dipole =
      (p3m.params.epsilon != P3M_EPSILON_METALLIC)
          ? std::make_optional(calc_dipole_moment(
                comm_cart, boost::combine(p_q_range, p_unfolded_pos_range)))
          : std::nullopt;
  auto const volume = box_geo.volume();
  auto const pref =
      4. * std::numbers::pi / volume / (2. * p3m.params.epsilon + 1.);

  /* === k-space force calculation  === */
  if (force_flag) {
    /* i*k differentiation */
    auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
    auto const &mesh_stop = p3m.mesh.size;
    auto const &offset = p3m.mesh.start;
    auto const wavevector = Utils::Vector3<FloatType>((2. * std::numbers::pi) *
                                                      box_geo.length_inv());
    auto indices = Utils::Vector3i{};
    auto index = std::size_t(0u);

    /* compute electric field */
    // Eq. (3.49) @cite deserno00b
    for_each_3d(mesh_start, mesh_stop, indices, [&]() {
      auto const rho_hat =
          std::complex<FloatType>(p3m.mesh.rs_scalar[2u * index + 0u],
                                  p3m.mesh.rs_scalar[2u * index + 1u]);
      auto const phi_hat = p3m.g_force[index] * rho_hat;

      for (int d = 0; d < 3; d++) {
        /* direction in r-space: */
        int d_rs = (d + p3m.mesh.ks_pnum) % 3;
        /* directions */
        auto const k = FloatType(p3m.d_op[d_rs][indices[d] + offset[d]]) *
                       wavevector[d_rs];

        /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
        p3m.mesh.rs_fields[d_rs][2u * index + 0u] = -k * phi_hat.imag();
        p3m.mesh.rs_fields[d_rs][2u * index + 1u] = +k * phi_hat.real();
      }

      ++index;
    });

    auto const check_residuals =
        not p3m.params.tuning and check_complex_residuals;
    p3m.fft->check_complex_residuals = check_residuals;
    for (auto &rs_mesh : p3m.fft_buffers->get_vector_mesh()) {
      p3m.fft->backward_fft(rs_mesh);
    }
    p3m.fft_buffers->perform_vector_halo_spread();
    p3m.fft->check_complex_residuals = false;

    auto const force_prefac = prefactor / volume;
    Utils::integral_parameter<int, AssignForces, 1, 7>(
        p3m.params.cao, p3m, force_prefac,
        boost::combine(p_q_range, p_force_range));

    // add dipole forces
    // Eq. (3.19) @cite deserno00b
    if (box_dipole) {
      auto const dm = prefactor * pref * box_dipole.value();
      for (auto zipped : boost::combine(p_q_range, p_force_range)) {
        auto p_q = boost::get<0>(zipped);
        auto &p_force = boost::get<1>(zipped);
        p_force -= p_q * dm;
      }
    }
  }

  /* === k-space energy calculation  === */
  if (energy_flag or npt_flag) {
    auto node_energy = 0.;
    auto const mesh_length = Utils::product(p3m.mesh.size);
    for (int i = 0; i < mesh_length; i++) {
      // Use the energy optimized influence function for energy!
      // Eq. (3.40) @cite deserno00b
      node_energy +=
          double(p3m.g_energy[i] * (Utils::sqr(p3m.mesh.rs_scalar[2 * i + 0]) +
                                    Utils::sqr(p3m.mesh.rs_scalar[2 * i + 1])));
    }
    node_energy /= 2. * volume;

    auto energy = 0.;
    boost::mpi::reduce(comm_cart, node_energy, energy, std::plus<>(), 0);
    if (this_node == 0) {
      /* self energy correction */
      // Eq. (3.8) @cite deserno00b
      energy -= p3m.sum_q2 * p3m.params.alpha * std::numbers::inv_sqrtpi;
      /* net charge correction */
      // Eq. (3.11) @cite deserno00b
      energy -= p3m.square_sum_q * std::numbers::pi /
                (2. * volume * Utils::sqr(p3m.params.alpha));
      /* dipole correction */
      // Eq. (3.9) @cite deserno00b
      if (box_dipole) {
        energy += pref * box_dipole.value().norm2();
      }
    }
    energy *= prefactor;
#ifdef NPT
    if (npt_flag) {
      npt_add_virial_contribution(energy);
    }
#endif
    if (not energy_flag) {
      energy = 0.;
    }
    return energy;
  }

  return 0.;
}

template <typename FloatType, Arch Architecture>
class CoulombTuningAlgorithm : public TuningAlgorithm {
  p3m_data_struct_coulomb<FloatType> &p3m;
  double m_mesh_density_min = -1., m_mesh_density_max = -1.;
  // indicates if mesh should be tuned
  bool m_tune_mesh = false;

protected:
  P3MParameters &get_params() override { return p3m.params; }

public:
  CoulombTuningAlgorithm(System::System &system, auto &input_p3m,
                         double prefactor, int timings)
      : TuningAlgorithm(system, prefactor, timings), p3m{input_p3m} {}

  void on_solver_change() const override { m_system.on_coulomb_change(); }

  void setup_logger(bool verbose) override {
    auto const &box_geo = *m_system.box_geo;
#ifdef CUDA
    auto const on_gpu = Architecture == Arch::GPU;
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

  std::optional<std::string>
  layer_correction_veto_r_cut(double r_cut) const override {
    auto const &solver = m_system.coulomb.impl->solver;
    if (auto actor = get_actor_by_type<ElectrostaticLayerCorrection>(solver)) {
      return actor->veto_r_cut(r_cut);
    }
    return {};
  }

  std::tuple<double, double, double, double>
  calculate_accuracy(Utils::Vector3i const &mesh, int cao,
                     double r_cut_iL) const override {

    auto const &box_geo = *m_system.box_geo;
    double alpha_L, rs_err, ks_err;

    /* calc maximal real space error for setting */
    rs_err = p3m_real_space_error(m_prefactor, r_cut_iL, p3m.sum_qpart,
                                  p3m.sum_q2, 0., box_geo.length());

    if (std::numbers::sqrt2 * rs_err > p3m.params.accuracy) {
      /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
      alpha_L = sqrt(log(std::numbers::sqrt2 * rs_err / p3m.params.accuracy)) /
                r_cut_iL;
    } else {
      /* even alpha=0 is ok, however, we cannot choose it since it kills the
         k-space error formula.
         Anyways, this very likely NOT the optimal solution */
      alpha_L = 0.1;
    }

    /* calculate real-space and k-space error for this alpha_L */
    rs_err = p3m_real_space_error(m_prefactor, r_cut_iL, p3m.sum_qpart,
                                  p3m.sum_q2, alpha_L, box_geo.length());
#ifdef CUDA
    if constexpr (Architecture == Arch::GPU) {
      if (this_node == 0) {
        ks_err =
            p3m_k_space_error_gpu(m_prefactor, mesh.data(), cao, p3m.sum_qpart,
                                  p3m.sum_q2, alpha_L, box_geo.length().data());
      }
      boost::mpi::broadcast(comm_cart, ks_err, 0);
    } else
#endif
      ks_err = p3m_k_space_error(m_prefactor, mesh, cao, p3m.sum_qpart,
                                 p3m.sum_q2, alpha_L, box_geo.length());

    return {Utils::Vector2d{rs_err, ks_err}.norm(), rs_err, ks_err, alpha_L};
  }

  void determine_mesh_limits() override {
    auto const &box_geo = *m_system.box_geo;
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
        for (auto i : {1, 2}) {
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
    auto const &box_geo = *m_system.box_geo;
    auto const &solver = m_system.coulomb.impl->solver;
    auto tuned_params = TuningAlgorithm::Parameters{};
    auto time_best = time_sentinel;
    auto mesh_density = m_mesh_density_min;
    while (mesh_density <= m_mesh_density_max) {
      auto trial_params = TuningAlgorithm::Parameters{};
      if (m_tune_mesh) {
        for (auto i : {0, 1, 2}) {
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
        if (has_actor_of_type<CoulombP3M>(solver)) {
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

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::tune() {
  auto &system = get_system();
  auto const &box_geo = *system.box_geo;
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
      CoulombTuningAlgorithm<FloatType, Architecture> parameters(
          system, p3m, prefactor, tune_timings);
      parameters.setup_logger(tune_verbose);
      // parameter ranges
      parameters.determine_mesh_limits();
      parameters.determine_r_cut_limits();
      parameters.determine_cao_limits(7);
      // run tuning algorithm
      parameters.tune();
      m_is_tuned = true;
      system.on_coulomb_change();
    } catch (...) {
      p3m.params.tuning = false;
      throw;
    }
  }
  init();
}

void CoulombP3M::sanity_checks_boxl() const {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  for (auto i = 0u; i < 3u; i++) {
    /* check k-space cutoff */
    if (p3m_params.cao_cut[i] >= box_geo.length_half()[i]) {
      std::stringstream msg;
      msg << "P3M_init: k-space cutoff " << p3m_params.cao_cut[i]
          << " is larger than half of box dimension " << box_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
    if (p3m_params.cao_cut[i] >= local_geo.length()[i]) {
      std::stringstream msg;
      msg << "P3M_init: k-space cutoff " << p3m_params.cao_cut[i]
          << " is larger than local box dimension " << local_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
  }

  if (p3m_params.epsilon != P3M_EPSILON_METALLIC) {
    if ((box_geo.length()[0] != box_geo.length()[1]) or
        (box_geo.length()[1] != box_geo.length()[2]) or
        (p3m_params.mesh[0] != p3m_params.mesh[1]) or
        (p3m_params.mesh[1] != p3m_params.mesh[2])) {
      throw std::runtime_error(
          "CoulombP3M: non-metallic epsilon requires cubic box");
    }
  }
}

void CoulombP3M::sanity_checks_periodicity() const {
  auto const &box_geo = *get_system().box_geo;
  if (!box_geo.periodic(0) or !box_geo.periodic(1) or !box_geo.periodic(2)) {
    throw std::runtime_error(
        "CoulombP3M: requires periodicity (True, True, True)");
  }
}

void CoulombP3M::sanity_checks_cell_structure() const {
  auto const &local_geo = *get_system().local_geo;
  if (local_geo.cell_structure_type() != CellStructureType::REGULAR and
      local_geo.cell_structure_type() != CellStructureType::HYBRID) {
    throw std::runtime_error(
        "CoulombP3M: requires the regular or hybrid decomposition cell system");
  }
  if (::communicator.size > 1 and
      local_geo.cell_structure_type() == CellStructureType::HYBRID) {
    throw std::runtime_error(
        "CoulombP3M: does not work with the hybrid decomposition cell system, "
        "if using more than one MPI node");
  }
}

void CoulombP3M::sanity_checks_node_grid() const {
  auto const &node_grid = ::communicator.node_grid;
  if (node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    throw std::runtime_error(
        "CoulombP3M: node grid must be sorted, largest first");
  }
}

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::scaleby_box_l() {
  auto const &box_geo = *get_system().box_geo;
  p3m.params.r_cut = p3m.params.r_cut_iL * box_geo.length()[0];
  p3m.params.alpha = p3m.params.alpha_L * box_geo.length_inv()[0];
  p3m.params.recalc_a_ai_cao_cut(box_geo.length());
  p3m.local_mesh.recalc_ld_pos(p3m.params);
  sanity_checks_boxl();
  calc_influence_function_force();
  calc_influence_function_energy();
}

#ifdef CUDA
template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::add_long_range_forces_gpu(
    ParticleRange const &particles) {
  if constexpr (Architecture == Arch::GPU) {
#ifdef NPT
    if (get_system().propagation->integ_switch == INTEG_METHOD_NPT_ISO) {
      auto const energy = long_range_energy(particles);
      npt_add_virial_contribution(energy);
    }
#else
    static_cast<void>(particles);
#endif
    if (this_node == 0) {
      auto &gpu = get_system().gpu;
      p3m_gpu_add_farfield_force(*m_gpu_data, gpu, prefactor,
                                 gpu.n_particles());
    }
  }
}

/* Initialize the CPU kernels.
 * This operation is time-consuming and sets up data members
 * that are only relevant for ELC force corrections, since the
 * GPU implementation uses CPU kernels to compute energies.
 */
template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::init_gpu_kernels() {
  if constexpr (Architecture == Arch::GPU) {
    auto &system = get_system();
    if (has_actor_of_type<ElectrostaticLayerCorrection>(
            system.coulomb.impl->solver)) {
      init_cpu_kernels();
    }
    p3m_gpu_init(m_gpu_data, p3m.params.cao, p3m.params.mesh, p3m.params.alpha,
                 system.box_geo->length(), system.gpu.n_particles());
  }
}

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::request_gpu() const {
  if constexpr (Architecture == Arch::GPU) {
    auto &gpu_particle_data = get_system().gpu;
    gpu_particle_data.enable_property(GpuParticleData::prop::force);
    gpu_particle_data.enable_property(GpuParticleData::prop::q);
    gpu_particle_data.enable_property(GpuParticleData::prop::pos);
  }
}
#endif // CUDA

template struct CoulombP3MImpl<float, Arch::CPU>;
template struct CoulombP3MImpl<float, Arch::GPU>;
template struct CoulombP3MImpl<double, Arch::CPU>;
template struct CoulombP3MImpl<double, Arch::GPU>;

#endif // P3M
