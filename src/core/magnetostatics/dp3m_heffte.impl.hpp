/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_DP3M

#include "magnetostatics/dp3m.hpp"

#include "magnetostatics/dipoles.hpp"
#include "short_range_cabana.hpp"

#include "magnetostatics/dp3m_heffte.hpp" // must be included after dipoles.hpp

#include "fft/fft.hpp"
#include "p3m/P3MFFT.hpp"
#include "p3m/TuningAlgorithm.hpp"
#include "p3m/TuningLogger.hpp"
#include "p3m/common.hpp"
#include "p3m/field_layout_helpers.hpp"
#include "p3m/influence_function_dipolar.hpp"
#include "p3m/interpolation.hpp"
#include "p3m/math.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrators/Propagation.hpp"
#include "kokkos_helpers.hpp"
#include "npt.hpp"
#include "system/System.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iterator>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
#ifndef NDEBUG
template <typename T>
bool heffte_almost_equal(T const &value, T const &reference) {
  auto const diff = std::abs(value - reference);
  using FT = std::remove_cvref_t<decltype(diff)>;
  auto constexpr atol = std::is_same_v<FT, float> ? FT{2e-4} : FT{1e-6};
  auto constexpr rtol = std::is_same_v<FT, float> ? FT{5e-5} : FT{1e-5};
  auto const non_zero = std::abs(reference) != FT{0};
  return (diff < atol) or (non_zero and (diff / std::abs(reference) < rtol));
}
#endif // not NDEBUG
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture,
                      FFTConfig>::count_magnetic_particles() {
  auto local_n = std::size_t{0u};
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

inline double dp3m_k_space_error(double box_size, int mesh, int cao,
                                 std::size_t n_c_part, double sum_q2,
                                 double alpha_L);

inline double dp3m_real_space_error(double box_size, double r_cut_iL,
                                    std::size_t n_c_part, double sum_q2,
                                    double alpha_L);

/** Compute the value of alpha through a bisection method.
 *  Based on eq. (33) @cite wang01a.
 */
double dp3m_rtbisection(double box_size, double r_cut_iL, std::size_t n_c_part,
                        double sum_q2, double x1, double x2, double xacc,
                        double tuned_accuracy);

template <typename FloatType, Arch Architecture, class FFTConfig>
double DipolarP3MHeffte<FloatType, Architecture,
                        FFTConfig>::calc_average_self_energy_k_space() const {
  auto const &box_geo = *get_system().box_geo;
  auto const node_phi =
      grid_influence_function_self_energy<FloatType, P3M_BRILLOUIN>(
          dp3m.params, dp3m.mesh.start, dp3m.mesh.stop, dp3m.g_energy);

  double phi = 0.;
  boost::mpi::reduce(comm_cart, node_phi, phi, std::plus<>(), 0);
  phi /= 3. * box_geo.length()[0] *
         Utils::int_pow<3>(static_cast<double>(dp3m.params.mesh[0]));
  return phi * std::numbers::pi;
}

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::init_cpu_kernels() {
  assert(dp3m.params.mesh >= Utils::Vector3i::broadcast(1));
  assert(dp3m.params.cao >= p3m_min_cao and dp3m.params.cao <= p3m_max_cao);
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
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  dp3m.heffte.world_size = comm_cart.size();
  dp3m.heffte.fft =
      std::make_shared<P3MFFT<FloatType, Architecture, FFTConfig>>(
          nullptr, ::comm_cart, dp3m.params.mesh, dp3m.local_mesh.ld_no_halo,
          dp3m.local_mesh.ur_no_halo, ::communicator.node_grid);
  dp3m.resize_heffte_buffers();
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  dp3m.calc_differential_operator();

  /* fix box length dependent constants */
  scaleby_box_l();

  count_magnetic_particles();
}

namespace {
template <int cao> struct AssignDipole {
  void operator()(auto &dp3m, auto &cell_structure) {
    using DipolarP3MState = std::remove_reference_t<decltype(dp3m)>;
    using value_type = DipolarP3MState::value_type;
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    auto const &aosoa = cell_structure.get_aosoa();
    auto const &unique_particles = cell_structure.get_unique_particles();
    auto const n_part = cell_structure.count_local_particles();
    dp3m.inter_weights.zfill(n_part); // allocate buffer for parallel write
    kokkos_parallel_range_for<execution_space>(
        "InterpolateDipoles", std::size_t{0u}, n_part, [&](auto p_index) {
          auto constexpr memory_order = Utils::MemoryOrder::ROW_MAJOR;
          auto const tid = omp_get_thread_num();
          auto const p_pos = aosoa.get_span_at(aosoa.position, p_index);
          auto const dip = unique_particles.at(p_index)->calc_dip();
          auto const weights =
              p3m_calculate_interpolation_weights<cao, memory_order>(
                  p_pos, dp3m.params.ai, dp3m.local_mesh);
          dp3m.inter_weights.store_at(p_index, weights);
          p3m_interpolate(
              dp3m.local_mesh, weights, [&dip, tid, &dp3m](int ind, double w) {
                dp3m.rs_fields_kokkos(tid, 0u, ind) += value_type(w * dip[0u]);
                dp3m.rs_fields_kokkos(tid, 1u, ind) += value_type(w * dip[1u]);
                dp3m.rs_fields_kokkos(tid, 2u, ind) += value_type(w * dip[2u]);
              });
        });
    Kokkos::fence();
    int num_threads = execution_space().concurrency();
    kokkos_parallel_range_for<execution_space>(
        "ReduceInterpolatedDipoles", std::size_t{0}, dp3m.local_mesh.size,
        [&dp3m, num_threads](std::size_t const i) {
          for (int dir = 0; dir < 3; ++dir) {
            value_type acc{};
            for (int tid = 0; tid < num_threads; ++tid) {
              acc += dp3m.rs_fields_kokkos(tid, dir, i);
            }
            dp3m.mesh.rs_fields[dir][i] += acc;
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
            dp3m.heffte.rs_dipole_density[dir][i] += acc;
#endif
          }
        });
    Kokkos::fence();
  }
};
} // namespace

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::dipole_assign() {
  prepare_fft_mesh();

  Utils::integral_parameter<int, AssignDipole, p3m_min_cao, p3m_max_cao>(
      dp3m.params.cao, dp3m, *get_system().cell_structure);
}

namespace {
template <int cao> struct AssignTorques {
  void operator()(auto &dp3m, double pref, int d_rs,
                  CellStructure &cell_structure) const {

    assert(cao == dp3m.inter_weights.cao());
    using execution_space = Kokkos::DefaultHostExecutionSpace;

    auto const kernel = [d_rs, pref, &dp3m](auto const &dip,
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
                                            auto &p_dip_fld,
#endif
                                            auto &p_torque,
                                            std::size_t p_index) {
      auto const weights = dp3m.inter_weights.template load<cao>(p_index);
      Utils::Vector3d E{};
      p3m_interpolate(dp3m.local_mesh, weights,
                      [&E, &dp3m, d_rs](int ind, double w) {
                        // heFFTe data: dp3m.heffte.ks_scalar.real()
                        E[d_rs] += w * double(dp3m.mesh.rs_scalar[ind]);
                      });

      auto const torque = pref * vector_product(dip, E);
      auto access = p_torque.access();
      access(p_index, 0) -= torque[0];
      access(p_index, 1) -= torque[1];
      access(p_index, 2) -= torque[2];
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      auto const dip_fld = pref * E;
      auto access_dip_fld = p_dip_fld.access();
      access_dip_fld(p_index, 0) -= dip_fld[0];
      access_dip_fld(p_index, 1) -= dip_fld[1];
      access_dip_fld(p_index, 2) -= dip_fld[2];
#endif
    };

    auto const n_part = dp3m.inter_weights.size();
    auto const &unique_particles = cell_structure.get_unique_particles();
    cell_structure.mark_torque_replicas_dirty();
    auto scatter_torque = cell_structure.get_scatter_torque();
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    auto scatter_dip_fld = cell_structure.get_scatter_dip_fld();
#endif
    kokkos_parallel_range_for<execution_space>(
        "AssignTorques", std::size_t{0u}, n_part, [&](std::size_t p_index) {
          auto const &p = *unique_particles.at(p_index);
          if (p.dipm() != 0.) {
            kernel(p.calc_dip(),
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
                   scatter_dip_fld,
#endif
                   scatter_torque, p_index);
          }
        });
  }
};

template <int cao> struct AssignForcesDip {
  void operator()(auto &dp3m, double pref, int d_rs,
                  CellStructure &cell_structure) const {

    assert(cao == dp3m.inter_weights.cao());
    using execution_space = Kokkos::DefaultHostExecutionSpace;

    auto const kernel = [d_rs, pref, &dp3m](auto const &dip, auto &p_force,
                                            std::size_t p_index) {
      auto const weights = dp3m.inter_weights.template load<cao>(p_index);

      Utils::Vector3d E{};
      p3m_interpolate(dp3m.local_mesh, weights, [&E, &dp3m](int ind, double w) {
        // heFFTe data: dp3m.heffte.rs_B_fields
        E[0u] += w * double(dp3m.mesh.rs_fields[0u][ind]);
        E[1u] += w * double(dp3m.mesh.rs_fields[1u][ind]);
        E[2u] += w * double(dp3m.mesh.rs_fields[2u][ind]);
      });

      auto access = p_force.access();
      access(p_index, d_rs) += pref * (dip * E);
    };

    auto const n_part = dp3m.inter_weights.size();
    auto const &unique_particles = cell_structure.get_unique_particles();
    auto scatter_force = cell_structure.get_scatter_force();
    kokkos_parallel_range_for<execution_space>(
        "AssignForcesDip", std::size_t{0u}, n_part, [&](std::size_t p_index) {
          auto const &p = *unique_particles.at(p_index);
          if (p.dipm() != 0.) {
            kernel(p.calc_dip(), scatter_force, p_index);
          }
        });
  }
};
} // namespace

#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
template <typename FloatType, class FFTConfig>
void DipolarP3MState<FloatType, FFTConfig>::resize_heffte_buffers() {
  auto const rs_array_size =
      static_cast<std::size_t>(Utils::product(this->local_mesh.dim));
  auto const rs_array_size_no_halo =
      static_cast<std::size_t>(Utils::product(this->local_mesh.dim_no_halo));
  auto const fft_mesh_size =
      static_cast<std::size_t>(Utils::product(heffte.fft->ks_local_size()));
  for (auto d : {0u, 1u, 2u}) {
    heffte.rs_dipole_density[d].resize(rs_array_size);
    heffte.ks_dipole_density[d].resize(fft_mesh_size);
    heffte.rs_B_fields[d].resize(rs_array_size);
    heffte.rs_B_fields_no_halo[d].resize(rs_array_size_no_halo);
  }
  heffte.ks_B_field_storage.resize(fft_mesh_size);
  heffte.ks_scalar.resize(fft_mesh_size);
}
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS

/**
 * @brief Reciprocal-space virial for the dipolar Ewald/P3M sum. Obtained
 * via the same Nose-Klein strain-derivative method used for the Coulomb
 * case (@cite essmann95a eq. (2.7), \f$\Pi_{\textrm{rec}, \alpha, \beta}\f$),
 * applied to the dipolar structure
 * factor \f$Q(\vec k) = \vec M(\vec k)\cdot\vec k\f$ with
 * \f$\vec M(\vec k) = \sum_j \vec \mu_j \exp(i\vec k\cdot\vec r_j)\f$.
 * Unlike the charge structure factor, \f$Q(\vec k)\f$ depends on
 * \f$\vec k\f$ explicitly (not only through the phase factor), which
 * produces an extra cross term beyond the charge-case \f$k_a k_b\f$
 * envelope. This cross term is generally asymmetric in \f$(a,b)\f$: its
 * symmetric half, \f$k_a\Re[M_b Q^*] + k_b\Re[M_a Q^*]\f$, is the
 * dipole-dipole reciprocal-space pressure tensor eq. (46) in
 * @cite aguado03a (their \f$\vec h\f$, \f$\kappa\f$ correspond to
 * \f$\vec k\f$, \f$\alpha\f$ here), which only ever reports that
 * symmetrized form. The remaining antisymmetric half is not in that
 * reference -- it is the reciprocal-space image of the same
 * dipole-dipole torque that already makes the real-space virial
 * asymmetric (see @ref DipolarDirectSum::long_range_pressure and
 * @ref DipolarP3M::pair_force), derived here by differentiating
 * the reciprocal energy directly (via the strain parametrization
 * \f$H(\varepsilon)=LI+\varepsilon E_{ab}\f$) instead of presupposing a
 * symmetric result.
 *
 * Care is needed with the index convention: probing the strain component
 * \f$\varepsilon_{ab}\f$ (i.e. \f$H(\varepsilon)=LI+\varepsilon E_{ab}\f$)
 * yields \f$-\partial U/\partial\varepsilon_{ab} = r_b F_a\f$ for a pair
 * separation \f$\vec r\f$ and force \f$\vec F\f$ -- the *transpose* of the
 * \f$r_a F_b\f$ (@ref Utils::tensor_product "d (x) f") convention used by
 * the real-space term and by @ref DipolarDirectSum::long_range_pressure.
 * Concretely, differentiating \f$Q(\vec k)=\vec k\cdot\vec M\f$ gives a
 * cross-term contribution to \f$-\partial U/\partial\varepsilon_{ab}\f$
 * proportional to \f$k_a\Re[M_bQ^*]\f$; to match the \f$r_aF_b\f$
 * convention, this must be stored as the \f$(b,a)\f$ tensor component,
 * i.e. \f$\Pi_{ab}\f$ gets cross term \f$2k_b\Re[M_aQ^*]\f$ (indices
 * swapped relative to the strain probe that produced it).
 */
template <typename FloatType, Arch Architecture, class FFTConfig>
Utils::Vector9d
DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::long_range_pressure() {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const dipole_prefac = prefactor / Utils::product(dp3m.params.mesh);
  Utils::Vector9d node_k_space_pressure_tensor{};

  if (dp3m.sum_mu2 > 0.) {
    dipole_assign();
    dp3m.fft_buffers->perform_vector_halo_gather();
    for (auto &rs_mesh : dp3m.fft_buffers->get_vector_mesh()) {
      dp3m.fft->forward_fft(rs_mesh);
    }
    dp3m.update_mesh_views();

    auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
    auto local_index = Utils::Vector3i::broadcast(0);
    auto const wavevector = 2. * std::numbers::pi * box_geo.length_inv()[0];
    auto const half_alpha_inv_sq = Utils::sqr(1. / (2. * dp3m.params.alpha));

    auto index = std::size_t(0u);
    auto it_energy = dp3m.g_energy.begin();
    for_each_3d(mesh_start, dp3m.mesh.size, local_index, [&]() {
      auto constexpr KX = 2, KY = 0, KZ = 1;
      auto const shift = local_index + dp3m.mesh.start;
      auto const &d_op = dp3m.d_op[0u];
      auto const &mesh_dip = dp3m.mesh.rs_fields;
      auto const d_op_x = static_cast<FloatType>(d_op[shift[KX]]);
      auto const d_op_y = static_cast<FloatType>(d_op[shift[KY]]);
      auto const d_op_z = static_cast<FloatType>(d_op[shift[KZ]]);

      // Re(M(k)) and Re(Q(k)) = Re(M(k)).n, same as the energy kernel's `re`
      auto const Mx_re = mesh_dip[0u][index];
      auto const My_re = mesh_dip[1u][index];
      auto const Mz_re = mesh_dip[2u][index];
      auto const Q_re = Mx_re * d_op_x + My_re * d_op_y + Mz_re * d_op_z;
      ++index;
      // Im(M(k)) and Im(Q(k))
      auto const Mx_im = mesh_dip[0u][index];
      auto const My_im = mesh_dip[1u][index];
      auto const Mz_im = mesh_dip[2u][index];
      auto const Q_im = Mx_im * d_op_x + My_im * d_op_y + Mz_im * d_op_z;
      ++index;

      auto const nx = static_cast<double>(d_op[shift[KX]]);
      auto const ny = static_cast<double>(d_op[shift[KY]]);
      auto const nz = static_cast<double>(d_op[shift[KZ]]);
      auto const kx = nx * wavevector;
      auto const ky = ny * wavevector;
      auto const kz = nz * wavevector;
      auto const norm_sq = Utils::sqr(kx) + Utils::sqr(ky) + Utils::sqr(kz);
      if (norm_sq != 0.) {
        auto const g = static_cast<double>(*it_energy);
        auto const cell_energy =
            g * static_cast<double>(Utils::sqr(Q_re) + Utils::sqr(Q_im));
        auto const vterm = -2. * (1. / norm_sq + half_alpha_inv_sq);

        // g * Re(M_a(k) * Q(k)^*), a in {x, y, z}
        auto const Rx = g * static_cast<double>(Mx_re * Q_re + Mx_im * Q_im);
        auto const Ry = g * static_cast<double>(My_re * Q_re + My_im * Q_im);
        auto const Rz = g * static_cast<double>(Mz_re * Q_re + Mz_im * Q_im);

        // Full (generally asymmetric) tensor: Pi_ab = cell_energy * (delta_ab
        // + vterm * k_a * k_b) + 2 * k_b * R_a (note: indices of the cross
        // term are swapped relative to the strain probe that produces it --
        // see the class-level comment above for the derivation). The
        // symmetric combination (R_a k_b + R_b k_a)/2 recovers the
        // literature (Aguado & Madden, eq. 46) result; the leftover
        // antisymmetric part is the k-space image of the same
        // dipolar-torque signature that already makes the real-space term
        // asymmetric (see dipolar_direct_sum.cpp / dp3m.hpp).
        node_k_space_pressure_tensor[0u] +=
            cell_energy * (1. + vterm * kx * kx) + 2. * nx * Rx; /* xx */
        node_k_space_pressure_tensor[1u] +=
            cell_energy * vterm * kx * ky + 2. * ny * Rx; /* xy */
        node_k_space_pressure_tensor[2u] +=
            cell_energy * vterm * kx * kz + 2. * nz * Rx; /* xz */
        node_k_space_pressure_tensor[3u] +=
            cell_energy * vterm * ky * kx + 2. * nx * Ry; /* yx */
        node_k_space_pressure_tensor[4u] +=
            cell_energy * (1. + vterm * ky * ky) + 2. * ny * Ry; /* yy */
        node_k_space_pressure_tensor[5u] +=
            cell_energy * vterm * ky * kz + 2. * nz * Ry; /* yz */
        node_k_space_pressure_tensor[6u] +=
            cell_energy * vterm * kz * kx + 2. * nx * Rz; /* zx */
        node_k_space_pressure_tensor[7u] +=
            cell_energy * vterm * kz * ky + 2. * ny * Rz; /* zy */
        node_k_space_pressure_tensor[8u] +=
            cell_energy * (1. + vterm * kz * kz) + 2. * nz * Rz; /* zz */
      }
      std::advance(it_energy, 1);
    });
  }

  return node_k_space_pressure_tensor * dipole_prefac * std::numbers::pi *
         box_geo.length_inv()[0];
}

template <typename FloatType, Arch Architecture, class FFTConfig>
double DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::long_range_kernel(
    bool force_flag, bool energy_flag) {
  /* k-space energy */
  double energy = 0.;
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const dipole_prefac = prefactor / Utils::product(dp3m.params.mesh);

  auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
  auto local_index = Utils::Vector3i::broadcast(0);
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  auto constexpr r2c_dir = FFTConfig::r2c_dir;
  auto const rs_local_size = dp3m.heffte.fft->rs_local_size();
  auto const local_size = dp3m.heffte.fft->ks_local_size();
  auto local_size_full = local_size;
  if constexpr (FFTConfig::use_r2c) {
    local_size_full[r2c_dir] -= 1;
    local_size_full[r2c_dir] *= 2;
  }
  auto const local_origin = dp3m.heffte.fft->ks_local_ld_index();
#ifndef NDEBUG
  auto const line_stride = local_size_full[0];
  auto const plane_stride = local_size_full[0] * local_size_full[0];
#endif
  auto const &global_size = dp3m.params.mesh;
  auto const cutoff_left = 1 - local_origin[r2c_dir];
  auto const cutoff_right = global_size[r2c_dir] / 2 - local_origin[r2c_dir];
  auto &short_dim = local_index[r2c_dir];
  dp3m.resize_heffte_buffers();
#endif

  if (dp3m.sum_mu2 > 0.) {
    dipole_assign();
    dp3m.fft_buffers->perform_vector_halo_gather();
    for (auto &rs_mesh : dp3m.fft_buffers->get_vector_mesh()) {
      dp3m.fft->forward_fft(rs_mesh);
    }
    dp3m.update_mesh_views();

#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
    if (dp3m.heffte.world_size == 1) {
      // halo communication of real space dipoles density
      std::array<FloatType *, 3u> rs_fields = {
          {dp3m.heffte.rs_dipole_density[0u].data(),
           dp3m.heffte.rs_dipole_density[1u].data(),
           dp3m.heffte.rs_dipole_density[2u].data()}};
      dp3m.heffte.halo_comm.gather_grid(::comm_cart, rs_fields,
                                        dp3m.local_mesh.dim);

      for (auto dir : {0u, 1u, 2u}) {
        // get real-space dipoles density without ghost layers
        extract_block_into<Utils::MemoryOrder::ROW_MAJOR,
                           FFTConfig::r_space_order>(
            dp3m.rs_field_no_halo_kokkos.data(),
            dp3m.heffte.rs_dipole_density[dir], dp3m.local_mesh.dim,
            dp3m.local_mesh.n_halo_ld,
            dp3m.local_mesh.dim - dp3m.local_mesh.n_halo_ur);
        // re-order data in row-major
        std::size_t index_row_major = 0u;
        for_each_3d_order<FFTConfig::k_space_order>(
            mesh_start, rs_local_size, local_index, [&]() {
              auto constexpr KX = 1, KY = 2, KZ = 0;
              auto const index = local_index[KZ] +
                                 rs_local_size[0] * local_index[KY] +
                                 Utils::sqr(rs_local_size[0]) * local_index[KX];
              dp3m.rs_field_no_halo_reorder_kokkos(index_row_major) =
                  dp3m.rs_field_no_halo_kokkos(index);
              ++index_row_major;
            });
        dp3m.heffte.fft->forward(dp3m.rs_field_no_halo_reorder_kokkos.data(),
                                 dp3m.heffte.ks_dipole_density[dir].data());
#ifndef NDEBUG
        if (not dp3m.params.tuning) {
          std::size_t index_row_major_r2c = 0u;
          for_each_3d_order<FFTConfig::k_space_order>(
              mesh_start, local_size, local_index, [&]() {
                if (not FFTConfig::use_r2c or (short_dim <= cutoff_right)) {
                  auto constexpr KX = 2, KY = 0, KZ = 1;
                  auto const index_fft_legacy = local_index[KZ] +
                                                line_stride * local_index[KY] +
                                                plane_stride * local_index[KX];
                  auto const old_value = std::complex<FloatType>{
                      dp3m.mesh.rs_fields[dir][2 * index_fft_legacy],
                      dp3m.mesh.rs_fields[dir][2 * index_fft_legacy + 1]};
                  auto const &new_value =
                      dp3m.heffte.ks_dipole_density[dir][index_row_major_r2c];
                  assert(heffte_almost_equal(new_value, old_value));
                  ++index_row_major_r2c;
                }
              });
        }
#endif // not NDEBUG
      }
    }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  }

  /* === k-space energy calculation  === */
  if (energy_flag) {
    /*********************
       Dipolar energy
    **********************/
    if (dp3m.sum_mu2 > 0.) {
      /* i*k differentiation for dipolar gradients:
       * |(\Fourier{\vect{mu}}(k)\cdot \vect{k})|^2 */

      auto index = std::size_t(0u);
      auto it_energy = dp3m.g_energy.begin();
      auto node_energy = 0.;
      for_each_3d(mesh_start, dp3m.mesh.size, local_index, [&]() {
        auto constexpr KX = 2, KY = 0, KZ = 1;
        auto const shift = local_index + dp3m.mesh.start;
        auto const &d_op = dp3m.d_op[0u];
        auto const &mesh_dip = dp3m.mesh.rs_fields;
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
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
      if (dp3m.heffte.world_size == 1) {
        [[maybe_unused]] auto node_energy_heffte = 0.;
        std::size_t index_row_major_r2c = 0u;
        for_each_3d_order<FFTConfig::k_space_order>(
            mesh_start, local_size, local_index, [&]() {
              if (not FFTConfig::use_r2c or (short_dim <= cutoff_right)) {
                auto const global_index = local_origin + local_index;
                auto const &mesh_dip = dp3m.heffte.ks_dipole_density;
                auto const cell_field =
                    mesh_dip[0u][index_row_major_r2c] *
                        FloatType(dp3m.d_op[0u][global_index[1u]]) +
                    mesh_dip[1u][index_row_major_r2c] *
                        FloatType(dp3m.d_op[1u][global_index[2u]]) +
                    mesh_dip[2u][index_row_major_r2c] *
                        FloatType(dp3m.d_op[2u][global_index[0u]]);
                auto cell_energy = static_cast<double>(
                    dp3m.heffte.g_energy[index_row_major_r2c] *
                    std::norm(cell_field));
                if (FFTConfig::use_r2c and (short_dim >= cutoff_left and
                                            short_dim <= cutoff_right - 1)) {
                  // k-space symmetry: double counting except in the first and
                  // last planes of the short dimension; although the wavevector
                  // points in the opposite direction in the redundant region of
                  // k-space, the product of two components of the wavevector
                  // cancels out the negative sign
                  cell_energy *= 2.;
                }
                node_energy_heffte += cell_energy;
              }
              ++index_row_major_r2c;
            });
        assert(heffte_almost_equal(static_cast<FloatType>(node_energy_heffte),
                                   static_cast<FloatType>(node_energy)));
      }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
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
      dp3m.ks_scalar.resize(dp3m.local_mesh.size);
      /* fill in ks_scalar array for torque calculation */
      {
        auto index{std::size_t(0u)};
        auto it_energy = dp3m.g_energy.begin();
        auto it_ks_scalar = dp3m.ks_scalar.begin();
        for_each_3d(mesh_start, dp3m.mesh.size, local_index, [&]() mutable {
          auto constexpr KX = 2, KY = 0, KZ = 1;
          auto const shift = local_index + dp3m.mesh.start;
          auto const &d_op = dp3m.d_op[0u];
          auto const &mesh_dip = dp3m.mesh.rs_fields;
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
          *it_ks_scalar = *it_energy * std::complex<FloatType>{re, im};
          std::advance(it_energy, 1);
          std::advance(it_ks_scalar, 1);
        });
      }
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
      if (dp3m.heffte.world_size == 1) {
        std::size_t index_row_major_r2c = 0u;
        for_each_3d_order<FFTConfig::k_space_order>(
            mesh_start, local_size, local_index, [&]() {
              if (not FFTConfig::use_r2c or (short_dim <= cutoff_right)) {
                auto const global_index = local_origin + local_index;
                auto const &mesh_dip = dp3m.heffte.ks_dipole_density;
                dp3m.heffte.ks_scalar[index_row_major_r2c] =
                    dp3m.heffte.g_energy[index_row_major_r2c] *
                    (mesh_dip[0u][index_row_major_r2c] *
                         FloatType(dp3m.d_op[0u][global_index[1u]]) +
                     mesh_dip[1u][index_row_major_r2c] *
                         FloatType(dp3m.d_op[1u][global_index[2u]]) +
                     mesh_dip[2u][index_row_major_r2c] *
                         FloatType(dp3m.d_op[2u][global_index[0u]]));
#ifndef NDEBUG
                if (not dp3m.params.tuning) {
                  auto constexpr KX = 2, KY = 0, KZ = 1;
                  auto const index_fft_legacy = local_index[KZ] +
                                                line_stride * local_index[KY] +
                                                plane_stride * local_index[KX];
                  assert(heffte_almost_equal(
                      dp3m.heffte.ks_scalar[index_row_major_r2c],
                      dp3m.ks_scalar[index_fft_legacy]));
                }
#endif // not NDEBUG
                ++index_row_major_r2c;
              }
            });
      }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS

      /* Torque component loop */
      for (int d = 0; d < 3; d++) {
        auto it_ks_scalar = dp3m.ks_scalar.begin();
        auto index = 0u;
        for_each_3d(mesh_start, dp3m.mesh.size, local_index, [&]() {
          auto const &offset = dp3m.mesh.start;
          auto const &d_op = dp3m.d_op[0u];
          auto const d_op_val = FloatType(d_op[local_index[d] + offset[d]]);
          auto const &value = *it_ks_scalar;
          dp3m.mesh.rs_scalar[index] = d_op_val * value.real();
          ++index;
          dp3m.mesh.rs_scalar[index] = d_op_val * value.imag();
          ++index;
          std::advance(it_ks_scalar, 1);
        });
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
        if (dp3m.heffte.world_size == 1) {
          unsigned int constexpr d_ks[3] = {2u, 0u, 1u};
          std::size_t index_row_major_r2c = 0u;
          for_each_3d_order<FFTConfig::k_space_order>(
              mesh_start, local_size, local_index, [&]() {
                if (not FFTConfig::use_r2c or (short_dim <= cutoff_right)) {
                  auto const global_index = local_origin + local_index;
                  auto const d_op_val =
                      FloatType(dp3m.d_op[d][global_index[d_ks[d]]]);
                  dp3m.heffte.ks_B_field_storage[index_row_major_r2c] =
                      d_op_val * dp3m.heffte.ks_scalar[index_row_major_r2c];
#ifndef NDEBUG
                  if (not dp3m.params.tuning) {
                    auto constexpr KX = 2, KY = 0, KZ = 1;
                    auto const index_fft_legacy =
                        local_index[KZ] + line_stride * local_index[KY] +
                        plane_stride * local_index[KX];
                    auto const old_value = std::complex<FloatType>{
                        dp3m.mesh.rs_scalar[2 * index_fft_legacy],
                        dp3m.mesh.rs_scalar[2 * index_fft_legacy + 1]};
                    auto const &new_value =
                        dp3m.heffte.ks_B_field_storage[index_row_major_r2c];
                    assert(heffte_almost_equal(new_value, old_value));
                  }
#endif // not NDEBUG
                  ++index_row_major_r2c;
                }
              });
          dp3m.heffte.fft->backward(dp3m.heffte.ks_B_field_storage.data(),
                                    dp3m.heffte.rs_B_fields_no_halo[d].data());
          // pad zeros around the B-field in real space for ghost layers,
          // writing straight into the persistent halo-sized buffer
          pad_with_zeros_discard_imag_into<FFTConfig::r_space_order,
                                           Utils::MemoryOrder::ROW_MAJOR>(
              dp3m.heffte.rs_B_fields[d].data(),
              std::span(dp3m.heffte.rs_B_fields_no_halo[d]),
              dp3m.local_mesh.dim_no_halo, dp3m.local_mesh.n_halo_ld,
              dp3m.local_mesh.n_halo_ur);
          // communicate ghost layers of the B-field in real space
          dp3m.heffte.halo_comm.spread_grid(::comm_cart,
                                            dp3m.heffte.rs_B_fields[d].data(),
                                            dp3m.local_mesh.dim);
        }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
        dp3m.fft->backward_fft(dp3m.fft_buffers->get_scalar_mesh());
        // communicate ghost layers of the B-field in real space
        dp3m.fft_buffers->perform_scalar_halo_spread();
        // assign torque component from mesh to particle
        auto const d_rs = (d + dp3m.mesh.ks_pnum) % 3;
        Utils::integral_parameter<int, AssignTorques, p3m_min_cao, p3m_max_cao>(
            dp3m.params.cao, dp3m, dipole_prefac * wavenumber, d_rs,
            *system.cell_structure);
      }

      /***************************
         DIPOLAR FORCES (k-space)
      ****************************/
      // Compute forces after torques because the algorithm below overwrites the
      // grids dp3m.mesh.rs_fields !
      // Note: I'll do here 9 inverse FFTs. By symmetry, we can reduce this
      // number to 6 !
      /* fill in ks_scalar array for force calculation */
      {
        auto it_force = dp3m.g_force.begin();
        auto it_ks_scalar = dp3m.ks_scalar.begin();
        std::size_t index = 0u;
        for_each_3d(mesh_start, dp3m.mesh.size, local_index, [&]() {
          auto constexpr KX = 2, KY = 0, KZ = 1;
          auto const shift = local_index + dp3m.mesh.start;
          auto const &d_op = dp3m.d_op[0u];
          auto const &mesh_dip = dp3m.mesh.rs_fields;
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
          *it_ks_scalar = {*it_force * im, *it_force * (-re)};
          std::advance(it_force, 1);
          std::advance(it_ks_scalar, 1);
        });
      }

#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
      if (dp3m.heffte.world_size == 1) {
        std::size_t index_row_major_r2c = 0u;
        for_each_3d_order<FFTConfig::k_space_order>(
            mesh_start, local_size, local_index, [&]() {
              if (not FFTConfig::use_r2c or (short_dim <= cutoff_right)) {
                auto const global_index = local_origin + local_index;
                auto const &mesh_dip = dp3m.heffte.ks_dipole_density;
                auto const value =
                    dp3m.heffte.g_force[index_row_major_r2c] *
                    (mesh_dip[0u][index_row_major_r2c] *
                         FloatType(dp3m.d_op[0u][global_index[1u]]) +
                     mesh_dip[1u][index_row_major_r2c] *
                         FloatType(dp3m.d_op[1u][global_index[2u]]) +
                     mesh_dip[2u][index_row_major_r2c] *
                         FloatType(dp3m.d_op[2u][global_index[0u]]));
                dp3m.heffte.ks_scalar[index_row_major_r2c] = {value.imag(),
                                                              -value.real()};
#ifndef NDEBUG
                if (not dp3m.params.tuning) {
                  auto constexpr KX = 2, KY = 0, KZ = 1;
                  auto const index_fft_legacy = local_index[KZ] +
                                                line_stride * local_index[KY] +
                                                plane_stride * local_index[KX];
                  assert(heffte_almost_equal(
                      dp3m.heffte.ks_scalar[index_row_major_r2c],
                      dp3m.ks_scalar[index_fft_legacy]));
                }
#endif // not NDEBUG
                ++index_row_major_r2c;
              }
            });
      }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS

      /* Force component loop */
      for (int d = 0; d < 3; d++) {
        std::size_t index = 0u;
        auto it_ks_scalar = dp3m.ks_scalar.begin();
        for_each_3d(mesh_start, dp3m.mesh.size, local_index, [&]() {
          auto constexpr KX = 2, KY = 0, KZ = 1;
          auto const shift = local_index + dp3m.mesh.start;
          auto const &d_op = dp3m.d_op[0u];
          auto const &mesh_dip = dp3m.mesh.rs_fields;
          auto const d_op_val = FloatType(d_op[shift[d]]);
          auto const f = *it_ks_scalar * d_op_val;
          mesh_dip[0u][index] = FloatType(d_op[shift[KX]]) * f.real();
          mesh_dip[1u][index] = FloatType(d_op[shift[KY]]) * f.real();
          mesh_dip[2u][index] = FloatType(d_op[shift[KZ]]) * f.real();
          ++index;
          mesh_dip[0u][index] = FloatType(d_op[shift[KX]]) * f.imag();
          mesh_dip[1u][index] = FloatType(d_op[shift[KY]]) * f.imag();
          mesh_dip[2u][index] = FloatType(d_op[shift[KZ]]) * f.imag();
          ++index;
          std::advance(it_ks_scalar, 1);
        });

#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
        if (dp3m.heffte.world_size == 1) {
          std::size_t index_row_major_r2c = 0u;
          for_each_3d_order<FFTConfig::k_space_order>(
              mesh_start, local_size, local_index, [&]() {
                if (not FFTConfig::use_r2c or (short_dim <= cutoff_right)) {
                  auto constexpr KX = 1, KY = 2, KZ = 0;
                  auto const global_index = local_origin + local_index;
                  auto const remapped_index =
                      local_index[KZ] + local_index[KY] * local_size[KZ] +
                      local_index[KX] * local_size[KZ] * local_size[KY];
                  auto const d_op_val =
                      FloatType(dp3m.d_op[d][global_index[d]]);
                  auto &mesh_dip = dp3m.heffte.ks_dipole_density;
                  mesh_dip[0u][index_row_major_r2c] =
                      FloatType(dp3m.d_op[d][global_index[2u]]) * d_op_val *
                      dp3m.heffte.ks_scalar[remapped_index];
                  mesh_dip[1u][index_row_major_r2c] =
                      FloatType(dp3m.d_op[d][global_index[0u]]) * d_op_val *
                      dp3m.heffte.ks_scalar[remapped_index];
                  mesh_dip[2u][index_row_major_r2c] =
                      FloatType(dp3m.d_op[d][global_index[1u]]) * d_op_val *
                      dp3m.heffte.ks_scalar[remapped_index];
#ifndef NDEBUG
                  if (not FFTConfig::use_r2c and not dp3m.params.tuning) {
                    auto const index_fft_legacy = local_index[2] +
                                                  line_stride * local_index[1] +
                                                  plane_stride * local_index[0];
                    for (int j = 0; j < 3; ++j) {
                      auto const old_value = std::complex<FloatType>{
                          dp3m.mesh.rs_fields[j][2 * index_fft_legacy],
                          dp3m.mesh.rs_fields[j][2 * index_fft_legacy + 1]};
                      auto const &new_value = mesh_dip[j][index_row_major_r2c];
                      assert(heffte_almost_equal(new_value, old_value));
                    }
                  }
#endif // not NDEBUG
                  ++index_row_major_r2c;
                }
              });
          for (int dir = 0u; dir < 3u; ++dir) {
            dp3m.heffte.fft->backward(
                dp3m.heffte.ks_dipole_density[dir].data(),
                dp3m.heffte.rs_B_fields_no_halo[dir].data());
            // pad zeros around the B-field in real space for ghost layers,
            // writing straight into the persistent halo-sized buffer
            pad_with_zeros_discard_imag_into<FFTConfig::r_space_order,
                                             Utils::MemoryOrder::ROW_MAJOR>(
                dp3m.heffte.rs_B_fields[d].data(),
                std::span(dp3m.heffte.rs_B_fields_no_halo[dir]),
                dp3m.local_mesh.dim_no_halo, dp3m.local_mesh.n_halo_ld,
                dp3m.local_mesh.n_halo_ur);
          }
          // communicate ghost layers of the B-field in real space
          auto rs_fields =
              std::array<FloatType *, 3u>{{dp3m.heffte.rs_B_fields[0u].data(),
                                           dp3m.heffte.rs_B_fields[1u].data(),
                                           dp3m.heffte.rs_B_fields[2u].data()}};
          dp3m.heffte.halo_comm.spread_grid(::comm_cart, rs_fields,
                                            dp3m.local_mesh.dim);
        }
#endif // ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
        for (auto &rs_mesh : dp3m.fft_buffers->get_vector_mesh()) {
          dp3m.fft->backward_fft(rs_mesh);
        }
        // communicate ghost layers of the B-field in real space
        dp3m.fft_buffers->perform_vector_halo_spread();
        // assign force component from mesh to particle
        auto const d_rs = (d + dp3m.mesh.ks_pnum) % 3;
        Utils::integral_parameter<int, AssignForcesDip, p3m_min_cao,
                                  p3m_max_cao>(
            dp3m.params.cao, dp3m, dipole_prefac * Utils::sqr(wavenumber), d_rs,
            *system.cell_structure);
      }
    } /* if (dp3m.sum_mu2 > 0) */
  } /* if (force_flag) */

  if (dp3m.params.epsilon != P3M_EPSILON_METALLIC) {
    auto const surface_term = calc_surface_term(force_flag, energy_flag);
    if (this_node == 0) {
      energy += surface_term;
    }
  }
#ifdef ESPRESSO_NPT
  if (force_flag and system.has_npt_enabled()) {
    // reuse the validated reciprocal-space pressure tensor (same one used by
    // the pressure observable) instead of an energy-proxy: unlike Coulomb,
    // the dipolar structure factor is not simply homogeneous in k, so energy
    // is not a valid substitute for the virial trace here (see
    // long_range_pressure())
    auto const pressure_tensor = long_range_pressure();
    get_system().npt_add_virial_contribution(
        pressure_tensor[0u] + pressure_tensor[4u] + pressure_tensor[8u]);
  }
#endif
  if (not energy_flag) {
    energy = 0.;
  }

  return energy;
}

template <typename FloatType, Arch Architecture, class FFTConfig>
double DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::calc_surface_term(
    bool force_flag, bool energy_flag) {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const particles = system.cell_structure->local_particles();
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
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      p.dip_fld() -= pref * box_dip;
#endif
      ip++;
    }
  }

  return energy;
}

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture,
                      FFTConfig>::calc_influence_function_force() {
  dp3m.g_force = grid_influence_function_dipolar<FloatType, 3, P3M_BRILLOUIN,
                                                 FFTConfig::k_space_order>(
      dp3m.params, dp3m.mesh.start, dp3m.mesh.stop,
      get_system().box_geo->length_inv());
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  if (dp3m.heffte.world_size == 1) {
    dp3m.heffte.g_force =
        grid_influence_function_dipolar<FloatType, 3, P3M_BRILLOUIN,
                                        FFTConfig::k_space_order>(
            dp3m.params, dp3m.heffte.fft->ks_local_ld_index(),
            dp3m.heffte.fft->ks_local_ur_index(),
            get_system().box_geo->length_inv());
    if constexpr (FFTConfig::use_r2c) {
      influence_function_r2c<FFTConfig::r2c_dir>(
          dp3m.heffte.g_force, dp3m.params.mesh,
          dp3m.heffte.fft->ks_local_size(),
          dp3m.heffte.fft->ks_local_ld_index());
    }
  }
#endif
}

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture,
                      FFTConfig>::calc_influence_function_energy() {
  dp3m.g_energy = grid_influence_function_dipolar<FloatType, 2, P3M_BRILLOUIN,
                                                  FFTConfig::k_space_order>(
      dp3m.params, dp3m.mesh.start, dp3m.mesh.stop,
      get_system().box_geo->length_inv());
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  if (dp3m.heffte.world_size == 1) {
    dp3m.heffte.g_energy =
        grid_influence_function_dipolar<FloatType, 2, P3M_BRILLOUIN,
                                        FFTConfig::k_space_order>(
            dp3m.params, dp3m.heffte.fft->ks_local_ld_index(),
            dp3m.heffte.fft->ks_local_ur_index(),
            get_system().box_geo->length_inv());
    if constexpr (FFTConfig::use_r2c) {
      influence_function_r2c<FFTConfig::r2c_dir>(
          dp3m.heffte.g_energy, dp3m.params.mesh,
          dp3m.heffte.fft->ks_local_size(),
          dp3m.heffte.fft->ks_local_ld_index());
    }
  }
#endif
}

template <typename FloatType, Arch Architecture, class FFTConfig>
class DipolarTuningAlgorithm : public TuningAlgorithm {
  DipolarP3MState<FloatType, FFTConfig> &dp3m;
  int m_mesh_max = -1, m_mesh_min = -1;
  std::pair<std::optional<int>, std::optional<int>> m_tune_limits;

public:
  DipolarTuningAlgorithm(System::System &system, decltype(dp3m) &input_dp3m,
                         double prefactor, int timings,
                         decltype(m_tune_limits) tune_limits)
      : TuningAlgorithm(system, prefactor, timings), dp3m{input_dp3m},
        m_tune_limits{std::move(tune_limits)} {}

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
      if (m_tune_limits.first) {
        m_mesh_min = *m_tune_limits.first;
      }
      if (m_tune_limits.second) {
        m_mesh_max = *m_tune_limits.second;
      }
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

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::tune() {
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
      DipolarTuningAlgorithm<FloatType, Architecture, FFTConfig> parameters(
          system, dp3m, prefactor, tuning.timings, tuning.limits);
      parameters.setup_logger(tuning.verbose);
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
inline auto dp3m_tune_aliasing_sums(Utils::Vector3i const &shift, int mesh,
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
inline double dp3m_k_space_error(double box_size, int mesh, int cao,
                                 std::size_t n_c_part, double sum_q2,
                                 double alpha_L) {

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
          if (d > 0. and std::fabs(d / alias1) > round_error_prec)
            he_q += d;
        }
      },
      [&values, &mesh_i, cotangent_sum](unsigned dim, int n) {
        values[dim] = cotangent_sum(n, mesh_i);
      });

  return 8. * Utils::sqr(std::numbers::pi) / 3. * sum_q2 *
         sqrt(he_q / static_cast<double>(n_c_part)) /
         Utils::int_pow<4>(box_size);
}

/** Calculate the value of the errors for the REAL part of the force in terms
 *  of the splitting parameter alpha of Ewald. Based on eq. (33) @cite wang01a.
 *
 *  Please note that in this more refined approach we don't use
 *  eq. (37), but eq. (33) which maintains all the powers in alpha.
 */
inline double dp3m_real_space_error(double box_size, double r_cut_iL,
                                    std::size_t n_c_part, double sum_q2,
                                    double alpha_L) {
  auto constexpr exp_min = -708.4; // for IEEE-compatible double
  double d_error_f, d_cc, d_dc, d_con;

  auto const d_rcut = r_cut_iL * box_size;
  auto const d_rcut2 = Utils::sqr(d_rcut);
  auto const d_rcut4 = Utils::sqr(d_rcut2);

  auto const d_a2 = Utils::sqr(alpha_L) / Utils::sqr(box_size);
  auto const exponent = -d_a2 * d_rcut2;
  auto const exp_term = (exponent < exp_min) ? 0. : std::exp(exponent);
  auto const d_c = sum_q2 * exp_term;

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
double dp3m_rtbisection(double box_size, double r_cut_iL, std::size_t n_c_part,
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

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture, FFTConfig>::scaleby_box_l() {
  auto const &box_geo = *get_system().box_geo;
  dp3m.params.r_cut = dp3m.params.r_cut_iL * box_geo.length()[0];
  dp3m.params.alpha = dp3m.params.alpha_L * box_geo.length_inv()[0];
  dp3m.params.recalc_a_ai_cao_cut(box_geo.length());
  dp3m.local_mesh.recalc_ld_pos(dp3m.params);
  sanity_checks_boxl();
  calc_influence_function_force();
  calc_influence_function_energy();
  dp3m.energy_correction = 0.;
#ifdef ESPRESSO_DP3M_HEFFTE_CROSS_CHECKS
  if (dp3m.heffte.world_size == 1) {
    dp3m.heffte.halo_comm.resize(::comm_cart, dp3m.local_mesh);
  }
#endif
}

template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture,
                      FFTConfig>::calc_energy_correction() {
  auto const &box_geo = *get_system().box_geo;
  auto const Ukp3m = calc_average_self_energy_k_space() * box_geo.volume();
  auto const Ewald_volume = Utils::int_pow<3>(dp3m.params.alpha_L);
  auto const Eself = -2. * Ewald_volume * std::numbers::inv_sqrtpi / 3.;
  dp3m.energy_correction =
      -dp3m.sum_mu2 * (Ukp3m + Eself + 2. * std::numbers::pi / 3.);
}

#ifdef ESPRESSO_NPT
template <typename FloatType, Arch Architecture, class FFTConfig>
void DipolarP3MHeffte<FloatType, Architecture,
                      FFTConfig>::npt_add_virial_contribution(double virial)
    const {
  get_system().npt_add_virial_contribution(virial);
}
#endif // ESPRESSO_NPT

#endif // ESPRESSO_DP3M
