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

#include "electrostatics/coulomb.hpp"
#include "electrostatics/elc.hpp"
#ifdef CUDA
#include "electrostatics/p3m_gpu_cuda.cuh"
#include "electrostatics/p3m_gpu_error.hpp"
#endif // CUDA

#include "electrostatics/p3m.impl.hpp"

#include "p3m/P3MFFT.hpp"
#include "p3m/TuningAlgorithm.hpp"
#include "p3m/TuningLogger.hpp"
#include "p3m/field_layout_helpers.hpp"
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
#include "aosoa_pack.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cell_system/particle_enumeration.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrators/Propagation.hpp"
#include "npt.hpp"
#include "p3m/send_mesh.hpp"
#include "particle_reduction.hpp"
#include "short_range_cabana.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>
#include <utils/serialization/array.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>

#ifdef SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#include <omp.h>
#endif

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
#include <vector>

template <typename FloatType>
std::complex<FloatType>
multiply_complex_by_imaginary(std::complex<FloatType> const &z, FloatType k) {
  // Perform the multiplication manually: (re + i*imag) * (i*k)
  return std::complex<FloatType>(-z.imag() * k, z.real() * k);
}

template <typename FloatType>
std::complex<FloatType>
multiply_complex_by_real(std::complex<FloatType> const &z, FloatType k) {
  // Perform the multiplication manually: (re + i*imag) * k
  return std::complex<FloatType>(z.real() * k, z.imag() * k);
}

template <typename FloatType>
FloatType complex_norm2(std::complex<FloatType> const &z) {
  return Utils::sqr(z.real()) + Utils::sqr(z.imag());
}

static bool is_node_grid_compatible_with_mesh(Utils::Vector3i const &node_grid,
                                              Utils::Vector3i const &mesh) {
  return mesh[0u] % node_grid[0u] == 0 and mesh[1u] % node_grid[1u] == 0 and
         mesh[2u] % node_grid[2u] == 0;
}

static auto get_size_from_shape(Utils::Vector3i const &shape) {
  return std::accumulate(shape.cbegin(), shape.cend(), std::size_t{1u},
                         [](std::size_t const acc, int const value) {
                           return acc * static_cast<std::size_t>(value);
                         });
}

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::count_charged_particles() {
  struct Res {
    std::size_t local_n = std::size_t{0u};
    double local_q = 0.0;
    double local_q2 = 0.0;
  };
  Reduction::AddPartialResultKernel<Res> kernel = [](Res &acc, auto const &p) {
    if (p.q() != 0.0) {
      acc.local_n++;
      acc.local_q2 += Utils::sqr(p.q());
      acc.local_q += p.q();
    }
  };

  Reduction::ReductionOp<Res> reduce = [](Res &a, Res const &b) {
    a.local_n += b.local_n;
    a.local_q += b.local_q;
    a.local_q2 += b.local_q2;
  };
  auto res = reduce_over_local_particles(*(get_system().cell_structure), kernel,
                                         reduce);

  boost::mpi::all_reduce(comm_cart, res.local_n, p3m.sum_qpart, std::plus<>());
  boost::mpi::all_reduce(comm_cart, res.local_q2, p3m.sum_q2, std::plus<>());
  boost::mpi::all_reduce(comm_cart, res.local_q, p3m.square_sum_q,
                         std::plus<>());
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
  p3m.g_force = grid_influence_function<FloatType, 1, P3M_BRILLOUIN>(
      p3m.params, p3m.fft->ks_local_ld_index(), p3m.fft->ks_local_ur_index(),
      get_system().box_geo->length_inv());
}

/** Calculate the influence function optimized for the energy and the
 *  self energy correction.
 */
template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::calc_influence_function_energy() {
  p3m.g_energy = grid_influence_function<FloatType, 0, P3M_BRILLOUIN>(
      p3m.params, p3m.fft->ks_local_ld_index(), p3m.fft->ks_local_ur_index(),
      get_system().box_geo->length_inv());
}

/** Aliasing sum used by @ref p3m_k_space_error. */
static auto p3m_tune_aliasing_sums(Utils::Vector3i const &shift,
                                   Utils::Vector3i const &mesh,
                                   Utils::Vector3d const &mesh_i, int cao,
                                   double alpha_L_i) {

  auto constexpr mesh_start = Utils::Vector3i::broadcast(-P3M_BRILLOUIN);
  auto constexpr mesh_stop = Utils::Vector3i::broadcast(P3M_BRILLOUIN + 1);
  auto constexpr exp_min = -708.4; // for IEEE-compatible double
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
        auto const exponent = -factor1 * norm_sq;
        auto const exp_limit = (exp_min + std::log(norm_sq)) / 2.;
        auto const ex = (exponent < exp_limit) ? 0. : std::exp(exponent);
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

  p3m.local_mesh.calc_local_ca_mesh(p3m.params, local_geo, skin, elc_layer);
  p3m.fft = std::make_shared<P3MFFT<FloatType>>(
      comm_cart, p3m.params.mesh, p3m.local_mesh.ld_no_halo,
      p3m.local_mesh.ur_no_halo, get_memory_layout());
  p3m.fft->set_preferred_kspace_decomposition(::communicator.node_grid);
  auto const rs_array_size = get_size_from_shape(p3m.local_mesh.dim);
  auto const rs_array_size_no_halo =
      get_size_from_shape(p3m.local_mesh.dim_no_halo);
  auto const fft_mesh_size = get_size_from_shape(p3m.fft->ks_local_size());
  p3m.rs_charge_density.resize(rs_array_size);
  p3m.ks_charge_density.resize(fft_mesh_size);
  for (auto d : {0u, 1u, 2u}) {
    p3m.rs_E_fields[d].resize(rs_array_size);
  }
  p3m.ks_E_fields_storage.resize(3u * fft_mesh_size);
  p3m.rs_E_fields_no_halo.resize(3u * rs_array_size_no_halo);
  p3m.calc_differential_operator();

  /* fix box length dependent constants */
  scaleby_box_l();

  count_charged_particles();
}

namespace {
template <int cao> struct AssignCharge {
  void operator()(auto &p3m, double q,
                  InterpolationWeights<cao> const &weights) {
    using value_type =
        typename std::remove_reference_t<decltype(p3m)>::value_type;
    p3m_interpolate(p3m.local_mesh, weights, [q, &p3m](int ind, double w) {
      p3m.rs_charge_density[ind] += value_type(w * q);
    });
  }

  void operator()(auto &p3m, double q, Utils::Vector3d const &real_pos,
                  p3m_interpolation_cache &inter_weights) {
    auto constexpr memory_order =
        std::remove_reference<decltype(p3m)>::type::memory_order;
    auto const weights = p3m_calculate_interpolation_weights<cao, memory_order>(
        real_pos, p3m.params.ai, p3m.local_mesh);
    inter_weights.store(weights);
    this->operator()(p3m, q, weights);
  }

  void operator()(auto &p3m, double q, Utils::Vector3d const &real_pos) {
    auto constexpr memory_order =
        std::remove_reference<decltype(p3m)>::type::memory_order;
    auto const weights = p3m_calculate_interpolation_weights<cao, memory_order>(
        real_pos, p3m.params.ai, p3m.local_mesh);
    this->operator()(p3m, q, weights);
  }

#ifdef SHARED_MEMORY_PARALLELISM
  void operator()(auto &p3m, auto &cell_structure) {
    auto const &aosoa = cell_structure.get_aosoa();
    auto const n_part = cell_structure.count_local_particles();
    p3m.inter_weights.zfill(n_part); // allocate buffer for parallel write
    kokkos_parallel_range_for(
        "InterpolationWeights", std::size_t{0u}, n_part, [&](auto p_index) {
          Utils::Vector3d const p_pos{aosoa.position(p_index, 0),
                                      aosoa.position(p_index, 1),
                                      aosoa.position(p_index, 2)};
          auto constexpr memory_order =
              std::remove_reference<decltype(p3m)>::type::memory_order;
          auto const weights =
              p3m_calculate_interpolation_weights<cao, memory_order>(
                  p_pos, p3m.params.ai, p3m.local_mesh);
          p3m.inter_weights.store_at(p_index, weights);
        });
    // serial part of charge assignment
    for (std::size_t p_index{0}; p_index < n_part; ++p_index) {
      auto const weights = p3m.inter_weights.template load<cao>(p_index);
      this->operator()(p3m, aosoa.charge(p_index), weights);
    }
  }
#else  // SHARED_MEMORY_PARALLELISM
  void operator()(auto &p3m, auto const &p_q_pos_range) {
    for (auto zipped : p_q_pos_range) {
      auto const p_q = boost::get<0>(zipped);
      if (p_q != 0.0) {
        auto const &p_pos = boost::get<1>(zipped);
        this->operator()(p3m, p_q, p_pos, p3m.inter_weights);
      }
    }
  }
#endif // SHARED_MEMORY_PARALLELISM
};
} // namespace

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::charge_assign(
    ParticleRange const &particles) {
  prepare_fft_mesh(true);

#ifdef SHARED_MEMORY_PARALLELISM
  Utils::integral_parameter<int, AssignCharge, 1, 7>(
      p3m.params.cao, p3m, *get_system().cell_structure);
#else  // SHARED_MEMORY_PARALLELISM
  auto p_q_range = ParticlePropertyRange::charge_range(particles);
  auto p_pos_range = ParticlePropertyRange::pos_range(particles);

  Utils::integral_parameter<int, AssignCharge, 1, 7>(
      p3m.params.cao, p3m, boost::combine(p_q_range, p_pos_range));
#endif // SHARED_MEMORY_PARALLELISM
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

namespace {
template <int cao> struct AssignForces {
  void operator()(auto &p3m, auto force_prefac,
#ifdef SHARED_MEMORY_PARALLELISM
                  CellStructure &cell_structure
#else
                  auto const &p_q_force_range
#endif
  ) const {

    assert(cao == p3m.inter_weights.cao());

    auto const kernel = [&p3m](auto pref, auto &p_force, std::size_t p_index) {
      auto const weights = p3m.inter_weights.template load<cao>(p_index);

      Utils::Vector3d force{};
      p3m_interpolate(p3m.local_mesh, weights,
                      [&force, &p3m](int ind, double w) {
                        force[0u] += w * double(p3m.rs_E_fields[0u][ind]);
                        force[1u] += w * double(p3m.rs_E_fields[1u][ind]);
                        force[2u] += w * double(p3m.rs_E_fields[2u][ind]);
                      });

#ifdef SHARED_MEMORY_PARALLELISM
      auto const thread_id = omp_get_thread_num();
      p_force(p_index, thread_id, 0) -= pref * force[0];
      p_force(p_index, thread_id, 1) -= pref * force[1];
      p_force(p_index, thread_id, 2) -= pref * force[2];
#else
      p_force -= pref * force;
#endif
    };

#ifdef SHARED_MEMORY_PARALLELISM
    auto const n_part = p3m.inter_weights.size();
    auto const &aosoa = cell_structure.get_aosoa();
    auto &local_force = cell_structure.get_local_force();
    kokkos_parallel_range_for(
        "AssignForces", std::size_t{0u}, n_part, [&](std::size_t p_index) {
          if (auto const pref = aosoa.charge(p_index) * force_prefac) {
            kernel(pref, local_force, p_index);
          }
        });
#else  // SHARED_MEMORY_PARALLELISM
    /* charged particle counter */
    std::size_t p_index{0ul};

    for (auto zipped : p_q_force_range) {
      auto p_q = boost::get<0>(zipped);
      if (p_q != 0.) {
        auto &p_force = boost::get<1>(zipped);
        kernel(p_q * force_prefac, p_force, p_index);
        ++p_index;
      }
    }
#endif // SHARED_MEMORY_PARALLELISM
  }
};
} // namespace

#ifdef SHARED_MEMORY_PARALLELISM
static auto calc_dipole_moment(boost::mpi::communicator const &comm,
                               auto const &cs, auto const &box_geo) {
  auto const local_dip = reduce_over_local_particles<Utils::Vector3d>(
      cs,
      [&box_geo](Utils::Vector3d &acc, Particle const &p) {
        acc += p.q() * box_geo.unfolded_position(p.pos(), p.image_box());
      },
      [](Utils::Vector3d &a, Utils::Vector3d const &b) { a = a + b; });
  return boost::mpi::all_reduce(comm, local_dip, std::plus<>());
}
#else  // SHARED_MEMORY_PARALLELISM
static auto calc_dipole_moment(boost::mpi::communicator const &comm,
                               auto const &p_q_unfolded_pos_range) {
  auto const local_dip =
      boost::accumulate(p_q_unfolded_pos_range, Utils::Vector3d{},
                        [](Utils::Vector3d const &dip, auto const &q_pos) {
                          auto const p_q = boost::get<0>(q_pos);
                          auto const &p_unfolded_pos = boost::get<1>(q_pos);
                          return dip + p_q * p_unfolded_pos;
                        });
  return boost::mpi::all_reduce(comm, local_dip, std::plus<>());
}
#endif // SHARED_MEMORY_PARALLELISM

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::kernel_ks_charge_density() {
  // halo communication of real space charge density
  p3m.halo_comm.gather_grid(comm_cart, p3m.rs_charge_density.data(),
                            p3m.local_mesh.dim);

  // get real space charge density without ghost layers
  auto charge_density_no_halos =
      extract_block<Utils::MemoryOrder::ROW_MAJOR,
                    Utils::MemoryOrder::COLUMN_MAJOR>(
          p3m.rs_charge_density, p3m.local_mesh.dim, p3m.local_mesh.n_halo_ld,
          p3m.local_mesh.dim - p3m.local_mesh.n_halo_ur);

  // Set up the FFT using the Heffte library.
  // This is in global mesh coordinates without any ghost layers
  // The memory layout has to be specified, so the parts of
  // the mesh held by each MPI rank are assembled correctly.
  p3m.fft->forward(charge_density_no_halos.data(),
                   p3m.ks_charge_density.data());
}

template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::kernel_rs_electric_field() {
  auto const mesh_start = p3m.fft->ks_local_ld_index();
  auto const mesh_stop = mesh_start + p3m.fft->ks_local_size();
  auto const &box_geo = *get_system().box_geo;

  // hold electric field in k-space
  std::array<std::span<std::complex<FloatType>>, 3> ks_E_fields;
  auto const fft_mesh_length = get_size_from_shape(mesh_stop - mesh_start);
  for (auto d : {0u, 1u, 2u}) {
    auto const offset = d * fft_mesh_length;
    auto const begin = p3m.ks_E_fields_storage.begin() + offset;
    ks_E_fields[d] = {begin, fft_mesh_length};
  }

  // i*k differentiation
  auto const wavevector =
      Utils::Vector3<FloatType>((2. * std::numbers::pi) * box_geo.length_inv());

  // compute electric field, Eq. (3.49) @cite deserno00b
  for_each_3d_lin<Utils::MemoryOrder::COLUMN_MAJOR>(
      mesh_start, mesh_stop,
      [&](const Utils::Vector3i &indices, int local_index) {
#ifdef ADDITIONAL_CHECKS
        assert(local_index ==
               Utils::get_linear_index<Utils::MemoryOrder::COLUMN_MAJOR>(
                   indices - mesh_start, p3m.fft->ks_local_size()));
#endif
        auto const phi_hat = multiply_complex_by_real(
            p3m.ks_charge_density[local_index], p3m.g_force[local_index]);

        for (auto d : {0u, 1u, 2u}) {
          // wave vector of the current mesh point
          auto const k = FloatType(p3m.d_op[d][indices[d]]) * wavevector[d];
          // electric field in k-space
          ks_E_fields[d][local_index] =
              multiply_complex_by_imaginary(phi_hat, k);
        }
      });

  // back-transform the k-space electric field to real space
  auto const size = p3m.local_mesh.ur_no_halo - p3m.local_mesh.ld_no_halo;
  auto const rs_mesh_size_no_halo = Utils::product(size);
  for (auto d : {0u, 1u, 2u}) {
    auto k_space = ks_E_fields[d].data();
    auto real_space = p3m.rs_E_fields_no_halo.data() + d * rs_mesh_size_no_halo;
    p3m.fft->backward(k_space, real_space);

    // add zeros around the E-field in real space to make room for ghost layers
    auto const offset = d * rs_mesh_size_no_halo;
    auto const begin = p3m.rs_E_fields_no_halo.begin() + offset;
    auto f = std::span<std::complex<FloatType>>(begin, rs_mesh_size_no_halo);
    p3m.rs_E_fields[d] =
        pad_with_zeros_discard_imag<Utils::MemoryOrder::COLUMN_MAJOR,
                                    Utils::MemoryOrder::ROW_MAJOR>(
            std::span<std::complex<FloatType>>(f), p3m.local_mesh.dim_no_halo,
            p3m.local_mesh.n_halo_ld, p3m.local_mesh.n_halo_ur);
  }

  // ghost communicate the boundary layers of the E-field in real space
  auto field_pointers = std::vector<FloatType *>{p3m.rs_E_fields[0u].data(),
                                                 p3m.rs_E_fields[1u].data(),
                                                 p3m.rs_E_fields[2u].data()};
  p3m.halo_comm.spread_grid(comm_cart, field_pointers, p3m.local_mesh.dim);
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
    kernel_ks_charge_density();

    auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
    auto const &mesh_stop = p3m.fft->ks_local_size();
    auto const &k_size = p3m.fft->ks_local_size();
    auto const half_alpha_inv_sq = Utils::sqr(1. / 2. / p3m.params.alpha);
    auto const wavevector = (2. * std::numbers::pi) * box_geo.length_inv();
    auto indices = Utils::Vector3i{};
    auto diagonal = 0.;

    for_each_3d(mesh_start, mesh_stop, indices, [&]() {
      auto const global_index = indices + p3m.fft->ks_local_ld_index();
      auto const kx = p3m.d_op[0u][global_index[0u]] * wavevector[0u];
      auto const ky = p3m.d_op[1u][global_index[1u]] * wavevector[1u];
      auto const kz = p3m.d_op[2u][global_index[2u]] * wavevector[2u];
      auto const norm_sq = Utils::sqr(kx) + Utils::sqr(ky) + Utils::sqr(kz);

      if (norm_sq != 0.) {
        auto const index = Utils::get_linear_index(
            indices, k_size, Utils::MemoryOrder::COLUMN_MAJOR);
        auto const node_k_space_energy = static_cast<double>(
            p3m.g_energy[index] * complex_norm2(p3m.ks_charge_density[index]));
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
  auto const npt_flag = force_flag and system.has_npt_enabled();
#else
  auto constexpr npt_flag = false;
#endif
  if (p3m.sum_qpart == 0u) {
    return 0.;
  }

  if (not has_actor_of_type<ElectrostaticLayerCorrection>(
          system.coulomb.impl->solver)) {
    charge_assign(particles);
  }

  kernel_ks_charge_density();

#ifdef SHARED_MEMORY_PARALLELISM
  auto &cell_structure = *system.cell_structure;
  auto const &local_force = cell_structure.get_local_force();
  auto const &aosoa = cell_structure.get_aosoa();
#else
  auto p_q_range = ParticlePropertyRange::charge_range(particles);
  auto p_force_range = ParticlePropertyRange::force_range(particles);
  auto p_unfolded_pos_range =
      ParticlePropertyRange::unfolded_pos_range(particles, box_geo);
#endif // SHARED_MEMORY_PARALLELISM

  // The dipole moment is only needed if we don't have metallic boundaries
  auto const box_dipole = (p3m.params.epsilon != P3M_EPSILON_METALLIC)
                              ? std::make_optional(calc_dipole_moment(
#ifdef SHARED_MEMORY_PARALLELISM
                                    comm_cart, cell_structure, box_geo))
#else
                                    comm_cart,
                                    boost::combine(p_q_range,
                                                   p_unfolded_pos_range)))
#endif
                              : std::nullopt;
  auto const volume = box_geo.volume();
  auto const pref =
      4. * std::numbers::pi / volume / (2. * p3m.params.epsilon + 1.);
  auto energy = 0.;

  /* === k-space force calculation  === */
  if (force_flag) {
    kernel_rs_electric_field();

    // assign particle forces
    auto const force_prefac = prefactor / volume;
#ifdef SHARED_MEMORY_PARALLELISM
    auto &particle_data = cell_structure;
#else
    auto const particle_data = boost::combine(p_q_range, p_force_range);
#endif
    Utils::integral_parameter<int, AssignForces, 1, 7>(
        p3m.params.cao, p3m, force_prefac, particle_data);

    // add dipole forces
    // Eq. (3.19) @cite deserno00b
    if (box_dipole) {
      auto const dm = prefactor * pref * box_dipole.value();
#ifdef SHARED_MEMORY_PARALLELISM
      auto const n_part = cell_structure.count_local_particles();
      kokkos_parallel_range_for(
          "AssignForcesBoxDipole", std::size_t{0u}, n_part,
          [&aosoa, &local_force, dm](auto p_index) {
            auto const thread_id = omp_get_thread_num();
            auto const q = aosoa.charge(p_index);
            local_force(p_index, thread_id, 0) -= q * dm[0];
            local_force(p_index, thread_id, 1) -= q * dm[1];
            local_force(p_index, thread_id, 2) -= q * dm[2];
          });
#else  // SHARED_MEMORY_PARALLELISM
      for (auto zipped : boost::combine(p_q_range, p_force_range)) {
        auto p_q = boost::get<0>(zipped);
        auto &p_force = boost::get<1>(zipped);
        p_force -= p_q * dm;
      }
#endif // SHARED_MEMORY_PARALLELISM
    }
  }

  /* === k-space energy calculation  === */
  if (energy_flag or npt_flag) {
    auto constexpr mesh_start = Utils::Vector3i::broadcast(0);
    auto const &mesh_stop = p3m.fft->ks_local_size();
    auto indices = Utils::Vector3i{};
    auto node_energy = 0.;
    for_each_3d(mesh_start, mesh_stop, indices, [&]() {
      auto const index = Utils::get_linear_index(
          indices, mesh_stop, Utils::MemoryOrder::COLUMN_MAJOR);
      // Use the energy optimized influence function for energy!
      // Eq. (3.40) @cite deserno00b
      node_energy += static_cast<double>(
          p3m.g_energy[index] * complex_norm2(p3m.ks_charge_density[index]));
    });
    node_energy /= 2. * volume;

    // add up energy contributions from all mpi ranks
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
      get_system().npt_add_virial_contribution(energy);
    }
#endif
    if (not energy_flag) {
      energy = 0.;
    }
  }

  return energy;
}

template <typename FloatType, Arch Architecture>
class CoulombTuningAlgorithm : public TuningAlgorithm {
  CoulombP3MState<FloatType> &p3m;
  double m_mesh_density_min = -1., m_mesh_density_max = -1.;
  // indicates if mesh should be tuned
  bool m_tune_mesh = false;
  std::pair<std::optional<int>, std::optional<int>> m_tune_limits;

protected:
  P3MParameters &get_params() override { return p3m.params; }

public:
  CoulombTuningAlgorithm(System::System &system, auto &input_p3m,
                         double prefactor, int timings,
                         decltype(m_tune_limits) tune_limits)
      : TuningAlgorithm(system, prefactor, timings), p3m{input_p3m},
        m_tune_limits{std::move(tune_limits)} {}

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

  std::optional<std::string> fft_decomposition_veto(
      Utils::Vector3i const &mesh_size_r_space) const override {
#ifdef CUDA
    if constexpr (Architecture == Arch::GPU) {
      return std::nullopt;
    }
#endif
    auto const [KX, KY, KZ] = p3m.fft->get_memory_layout();
    auto valid_decomposition = false;
    Utils::Vector3i mesh_size_k_space = {};
    boost::mpi::reduce(
        ::comm_cart, p3m.fft->ks_local_ur_index(), mesh_size_k_space,
        [](Utils::Vector3i const &lhs, Utils::Vector3i const &rhs) {
          return Utils::Vector3i{{std::max(lhs[0u], rhs[0u]),
                                  std::max(lhs[1u], rhs[1u]),
                                  std::max(lhs[2u], rhs[2u])}};
        },
        0);
    if (::this_node == 0) {
      auto const &node_grid = ::communicator.node_grid;
      valid_decomposition =
          (mesh_size_r_space[0u] == mesh_size_k_space[KX] and
           mesh_size_r_space[1u] == mesh_size_k_space[KY] and
           mesh_size_r_space[2u] == mesh_size_k_space[KZ] and
           is_node_grid_compatible_with_mesh(node_grid, mesh_size_r_space));
    }
    boost::mpi::broadcast(::comm_cart, valid_decomposition, 0);
    std::optional<std::string> retval{"conflict with FFT domain decomposition"};
    if (valid_decomposition) {
      retval = std::nullopt;
    }
    return retval;
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
      if (m_tune_limits.first or m_tune_limits.second) {
        auto const &box_l = box_geo.length();
        auto const dim = std::max({box_l[0], box_l[1], box_l[2]});
        if (m_tune_limits.first) {
          m_mesh_density_min = static_cast<double>(*m_tune_limits.first) / dim;
        }
        if (m_tune_limits.second) {
          m_mesh_density_max = static_cast<double>(*m_tune_limits.second) / dim;
        }
      }
      m_tune_mesh = true;
    } else {
      m_mesh_density_min = m_mesh_density_max = mesh_density;
      assert(p3m.params.mesh[0] >= 1);
      if (p3m.params.mesh[1] == -1 and p3m.params.mesh[2] == -1) {
        // determine the two missing values by rescaling by the box length
        for (auto i : {1u, 2u}) {
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
    auto current_mesh = p3m.params.mesh;
    if (m_tune_mesh) {
      for (auto i : {0u, 1u, 2u}) {
        current_mesh[i] =
            static_cast<int>(std::round(box_geo.length()[i] * mesh_density));
        // make the mesh even in all directions
        current_mesh[i] += current_mesh[i] % 2;
      }
    }

    while (mesh_density <= m_mesh_density_max) {
      auto trial_params = TuningAlgorithm::Parameters{};
      trial_params.mesh = current_mesh;
      trial_params.cao = cao_best;
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
      if (m_tune_mesh) {
        current_mesh += Utils::Vector3i::broadcast(2);
        mesh_density = current_mesh[0] / box_geo.length()[0];
      } else {
        break;
      }
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
          system, p3m, prefactor, tune_timings, tune_limits);
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
  p3m.halo_comm.resize(::comm_cart, p3m.local_mesh);
}

#ifdef CUDA
template <typename FloatType, Arch Architecture>
void CoulombP3MImpl<FloatType, Architecture>::add_long_range_forces_gpu(
    ParticleRange const &particles) {
  if constexpr (Architecture == Arch::GPU) {
#ifdef NPT
    if (get_system().has_npt_enabled()) {
      get_system().npt_add_virial_contribution(long_range_energy(particles));
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
