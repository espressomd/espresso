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

#ifdef ESPRESSO_DIPOLES

#include "magnetostatics/dipolar_direct_sum.hpp"
#include "magnetostatics/dipolar_direct_sum_kernels.hpp"

#include "BoxGeometry.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/mpi/iall_gatherv.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <cassert>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

/**
 * @brief Position and dipole moment of one particle.
 */
struct PosMom {
  Utils::Vector3d pos;
  Utils::Vector3d m;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar & pos & m;
  }
};

static auto gather_particle_data(BoxGeometry const &box_geo,
                                 ParticleRange const &particles) {
  auto const &comm = ::comm_cart;
  std::vector<Particle *> local_particles;
  std::vector<PosMom> local_posmom;
  std::vector<PosMom> all_posmom;
  std::vector<boost::mpi::request> reqs;

  local_particles.reserve(particles.size());
  local_posmom.reserve(particles.size());

  for (auto &p : particles) {
    if (p.dipm() != 0.0) {
      local_particles.emplace_back(&p);
      local_posmom.emplace_back(
          PosMom{box_geo.folded_position(p.pos()), p.calc_dip()});
    }
  }

  auto const local_size = static_cast<int>(local_posmom.size());
  std::vector<int> all_sizes;
  boost::mpi::all_gather(comm, local_size, all_sizes);

  auto const offset =
      std::accumulate(all_sizes.begin(), all_sizes.begin() + comm.rank(), 0);
  auto const total_size =
      std::accumulate(all_sizes.begin() + comm.rank(), all_sizes.end(), offset);

  if (comm.size() > 1) {
    all_posmom.resize(total_size);
    reqs = Utils::Mpi::iall_gatherv(comm, local_posmom.data(), local_size,
                                    all_posmom.data(), all_sizes.data());
  } else {
    std::swap(all_posmom, local_posmom);
  }

  return std::make_tuple(std::move(local_particles), std::move(all_posmom),
                         std::move(reqs), offset);
}

static auto get_n_cut(BoxGeometry const &box_geo, int n_replicas) {
  return n_replicas * Utils::Vector3i{static_cast<int>(box_geo.periodic(0)),
                                      static_cast<int>(box_geo.periodic(1)),
                                      static_cast<int>(box_geo.periodic(2))};
}

/**
 * Real-space image shifts n x box_l inside the |ncut| sphere. Index 0 is the
 * primary (zero) shift so self-interaction loops start at index 1.
 */
static std::vector<Utils::Vector3d>
make_image_shifts(Utils::Vector3i const &ncut, Utils::Vector3d const &box_l) {
  auto const ncut2 = ncut.norm2();
  std::vector<Utils::Vector3d> shifts;
  shifts.push_back({0., 0., 0.});
  for (int nx = -ncut[0]; nx <= ncut[0]; ++nx)
    for (int ny = -ncut[1]; ny <= ncut[1]; ++ny)
      for (int nz = -ncut[2]; nz <= ncut[2]; ++nz) {
        if (nx == 0 && ny == 0 && nz == 0)
          continue;
        if (nx * nx + ny * ny + nz * nz <= ncut2)
          shifts.push_back({nx * box_l[0], ny * box_l[1], nz * box_l[2]});
      }
  return shifts;
}

/**
 * @brief Calculate and add the interaction forces/torques to the particles.
 *
 * This employs a parallel N-square loop over all particle pairs.
 * The computation the partitioned into several steps so that the
 * communication latency can be hidden behind some local computation:
 *
 * 1. The local particle positions and momenta are packed into
 *    one array.
 * 2. The asynchronous distribution of the local arrays to all
 *    ranks is started.
 * 3. The interaction for the local pairs is started, here every
 *    pair is visited only once, and the force is added to both particles.
 * 4. Wait for the data from the other nodes.
 * 5. Calculate the interaction with the rest of the particles. Here
 *    every pair is visited twice (not necessarily on the same rank)
 *    so that no reduction of the forces is needed.
 *
 * Logically this is equivalent to the potential calculation
 * in @ref DipolarDirectSum::long_range_energy, which calculates
 * a naive N-square sum, but has better performance and scaling.
 */
void DipolarDirectSum::add_long_range_forces_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  /* Number of image boxes considered */
  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();

  auto const prefactor_local = prefactor;

  /* Raw pointer to the gathered AoS data; the local slice is populated before
   * wait_all, so Phase A may read it, and it outlives all fences. Safe to
   * capture by value in the Kokkos [=] lambdas. */
  auto const *pm = all_posmom.data();

  using execution_space = Kokkos::DefaultExecutionSpace;
  using ForceView =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using ScatterForce =
      Kokkos::Experimental::ScatterView<double *[3], Kokkos::LayoutRight>;
  ForceView local_force("dds_force", n_local);
  ForceView local_torque("dds_torque", n_local);
  ScatterForce scatter_force(local_force);
  ScatterForce scatter_torque(local_torque);

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto *local_particles_ptr = local_particles.data();
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  /* Phase A: local pairs. Each i owns its own force/torque accumulation
   * (written directly, unique owner, no race); the Newton's-third-law
   * partner-j contributions go through the ScatterView with a per-lane
   * scatter. */
  Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
  policy_local(std::size_t{0}, n_local);
  policy_local.set_chunk_size(64);
  Kokkos::parallel_for(
      "dds_local_pairs", policy_local, [=](std::size_t const i) {
        auto const gi = offset + i;
        auto const &pos_i = pm[gi].pos;
        auto const &m_i = pm[gi].m;
        PairForce fi{};

        /* (a) self-images (shifts[1..], primary excluded) */
        for (std::size_t s = 1; s < n_shifts; ++s)
          fi += pair_force(shifts_ptr[s], m_i, m_i);

        auto force_access = scatter_force.access();
        auto torque_access = scatter_torque.access();

        /* (b) pairs with j in (gi, offset + n_local) */
        for (auto j = gi + 1; j < offset + n_local; ++j) {
          auto const &pos_j = pm[j].pos;
          auto const &m_j = pm[j].m;
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          auto const jl = j - offset;
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const rn = d0 + shifts_ptr[s];
            auto const pf = pair_force(rn, m_i, m_j);
            fi.f += pf.f;
            fi.torque += pf.torque;
            /* Conservation of angular momentum mandates that
             * 0 = t_i + r_ij x F_ij + t_j */
            auto const torque_j = vector_product(pf.f, rn) - pf.torque;
            for (int c = 0; c < 3; ++c) {
              force_access(jl, c) -= pf.f[c];
              torque_access(jl, c) += torque_j[c];
            }
          }
        }
        /* (d) write i's own total directly (unique owner, no race) */
        local_particles_ptr[i]->force() += prefactor_local * fi.f;
        local_particles_ptr[i]->torque() += prefactor_local * fi.torque;
      });
  Kokkos::fence();

  /* Wait for remote data; the remote slices of all_posmom are now populated. */
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  /* Phase B: remote pairs (red [0, offset) + black [offset + n_local,
   * n_total)), visit-twice, no scatter; accumulate only i. */
  Kokkos::RangePolicy<execution_space> policy_remote(std::size_t{0}, n_local);
  Kokkos::parallel_for(
      "dds_remote_pairs", policy_remote, [=](std::size_t const i) {
        auto const gi = offset + i;
        auto const &pos_i = pm[gi].pos;
        auto const &m_i = pm[gi].m;
        PairForce fi{};

        /* Two remote ranges: red [0, offset) and black [offset + n_local,
         * n_total). Visit-twice (each remote pair is visited once per owning
         * rank), so only i accumulates; no scatter. */
        std::size_t const ranges[2][2] = {{std::size_t{0}, offset},
                                          {offset + n_local, n_total}};
        for (auto const &range : ranges) {
          auto const range_begin = range[0];
          auto const range_end = range[1];
          for (auto j = range_begin; j < range_end; ++j) {
            auto const &pos_j = pm[j].pos;
            auto const &m_j = pm[j].m;
            auto const d0 = with_replicas ? (pos_i - pos_j)
                                          : box_geo.get_mi_vector(pos_i, pos_j);
            for (std::size_t s = 0; s < n_shifts; ++s) {
              auto const rn = d0 + shifts_ptr[s];
              auto const pf = pair_force(rn, m_i, m_j);
              fi.f += pf.f;
              fi.torque += pf.torque;
            }
          }
        }
        local_particles_ptr[i]->force() += prefactor_local * fi.f;
        local_particles_ptr[i]->torque() += prefactor_local * fi.torque;
      });
  Kokkos::fence();

  /* Reduce the Newton's-third-law contributions and add to particles. */
  Kokkos::Experimental::contribute(local_force, scatter_force);
  Kokkos::Experimental::contribute(local_torque, scatter_torque);
  Kokkos::RangePolicy<execution_space> policy_reduce(std::size_t{0}, n_local);
  Kokkos::parallel_for(
      "dds_reduction", policy_reduce, [=](std::size_t const i) {
        local_particles_ptr[i]->force() +=
            prefactor_local * Utils::Vector3d{local_force(i, 0),
                                              local_force(i, 1),
                                              local_force(i, 2)};
        local_particles_ptr[i]->torque() +=
            prefactor_local * Utils::Vector3d{local_torque(i, 0),
                                              local_torque(i, 1),
                                              local_torque(i, 2)};
      });
  Kokkos::fence();

#ifdef ESPRESSO_NPT
  // As with DipolarP3M, the energy is not a valid
  // substitute for the virial trace for dipole-dipole interactions, so the
  // pressure tensor (see long_range_pressure()) is reused instead.
  if (system.has_npt_enabled()) {
    auto const pressure_tensor = long_range_pressure();
    get_system().npt_add_virial_contribution(
        pressure_tensor[0u] + pressure_tensor[4u] + pressure_tensor[8u]);
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (not m_is_gpu) {
    dipole_field_at_part_cpu();
  }
#endif
}

/**
 * @brief Calculate the interaction potential.
 *
 * This employs a parallel N-square loop over all particle pairs.
 */
double DipolarDirectSum::long_range_energy_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  /* Number of image boxes considered */
  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();

  /* Raw pointer to the gathered AoS data; the local slice is populated before
   * wait_all, so Phase A may read it, and it outlives all fences. Safe to
   * capture by value in the Kokkos [=] lambdas. */
  auto const *pm = all_posmom.data();

  using execution_space = Kokkos::DefaultExecutionSpace;

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  /* Phase A: local-upper triangular sum over j in [gi, offset + n_local),
   * i.e. the self-image energy (shifts[1..], primary excluded) plus the pairs
   * with j in (gi, offset + n_local). Computed from the local slice of the
   * gathered data while the remote data is still in flight. */
  double uA = 0.;
  Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
  policy_local(std::size_t{0}, n_local);
  policy_local.set_chunk_size(64);
  Kokkos::parallel_reduce(
      "dds_energy_local", policy_local,
      [=](std::size_t const i, double &u_local) {
        auto const gi = offset + i;
        auto const &pos_i = pm[gi].pos;
        auto const &m_i = pm[gi].m;

        /* (a) self-images (shifts[1..], primary excluded) */
        for (std::size_t s = 1; s < n_shifts; ++s)
          u_local += pair_potential(shifts_ptr[s], m_i, m_i);

        /* (b) pairs with j in (gi, offset + n_local) */
        for (auto j = gi + 1; j < offset + n_local; ++j) {
          auto const &pos_j = pm[j].pos;
          auto const &m_j = pm[j].m;
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          for (std::size_t s = 0; s < n_shifts; ++s)
            u_local += pair_potential(d0 + shifts_ptr[s], m_i, m_j);
        }
      },
      uA);

  /* Wait for remote data, fill the black slice. The red range [0, offset) is
   * never summed by the energy kernel (each pair is counted once on the rank
   * owning its lower index). */
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  /* Phase B: remote-black sum over j in [offset + n_local, n_total). No self
   * term and no primary exclusion; the range is entirely remote. */
  double uB = 0.;
  Kokkos::RangePolicy<execution_space> policy_remote(std::size_t{0}, n_local);
  Kokkos::parallel_reduce(
      "dds_energy_remote", policy_remote,
      [=](std::size_t const i, double &u_local) {
        auto const gi = offset + i;
        auto const &pos_i = pm[gi].pos;
        auto const &m_i = pm[gi].m;
        /* sum over j in [offset + n_local, n_total) */
        for (auto j = offset + n_local; j < n_total; ++j) {
          auto const &pos_j = pm[j].pos;
          auto const &m_j = pm[j].m;
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          for (std::size_t s = 0; s < n_shifts; ++s)
            u_local += pair_potential(d0 + shifts_ptr[s], m_i, m_j);
        }
      },
      uB);

  return prefactor * (uA + uB);
}

/**
 * @brief Calculate the dipolar pressure tensor.
 *
 * This employs a parallel N-square loop over all particle pairs.
 * The dipole-dipole force is not central (unlike Coulomb), so the pair
 * contribution to the virial only accounts for the force, not the torque,
 * matching the convention used by @ref DipolarP3M.
 */
Utils::Vector9d DipolarDirectSum::long_range_pressure() const {
  if (m_is_gpu) {
    runtimeWarningMsg() << "Pressure calculation not implemented for "
                           "DipolarDirectSum on GPU.";
    return Utils::Vector9d{};
  }

  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  /* Number of image boxes considered */
  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();

  /* Raw pointer to the gathered AoS data; the local slice is populated before
   * wait_all, so Phase A may read it, and it outlives all fences. Safe to
   * capture by value in the Kokkos [=] lambdas. */
  auto const *pm = all_posmom.data();

  using execution_space = Kokkos::DefaultExecutionSpace;

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  /* The pair contribution to the virial is r_n (x) F(r_n), reduced as a flat
   * 9-component array. The dipole-dipole force is not central, so only the
   * force enters the virial (not the torque). */

  /* Phase A: local upper-triangular sum over j in [gi, offset + n_local) --
   * the self-image term (shifts[1..], primary excluded) plus the pairs with
   * j in (gi, offset + n_local). Computed while remote data is in flight. */
  Utils::Vector9d pA{};
  Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>
  policy_local(std::size_t{0}, n_local);
  policy_local.set_chunk_size(64);
  auto const join_op = [](Utils::Vector9d &acc, Utils::Vector9d const &val) {
    acc += val;
  };
  auto reducerA = Reduction::make_kokkos_reducer<Utils::Vector9d>(
      // NOLINTNEXTLINE(bugprone-exception-escape)
      [=](std::size_t const i, Utils::Vector9d &psum) noexcept {
        auto const gi = offset + i;
        auto const &pos_i = pm[gi].pos;
        auto const &m_i = pm[gi].m;

        /* (a) self-images (shifts[1..], primary excluded) */
        for (std::size_t s = 1; s < n_shifts; ++s) {
          auto const rn = shifts_ptr[s];
          psum += Utils::flatten(
              Utils::tensor_product(rn, pair_force(rn, m_i, m_i).f));
        }

        /* (b) pairs with j in (gi, offset + n_local) */
        for (auto j = gi + 1; j < offset + n_local; ++j) {
          auto const &pos_j = pm[j].pos;
          auto const &m_j = pm[j].m;
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const rn = d0 + shifts_ptr[s];
            psum += Utils::flatten(
                Utils::tensor_product(rn, pair_force(rn, m_i, m_j).f));
          }
        }
      },
      join_op);
  Kokkos::parallel_reduce("dds_pressure_local", policy_local, reducerA, pA);

  /* Wait for remote data, fill the black slice. The red range [0, offset) is
   * never summed (each pair is counted once on the rank owning its lower
   * index). */
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  /* Phase B: remote-black sum over j in [offset + n_local, n_total). No self
   * term and no primary exclusion -- the range is entirely remote. */
  Utils::Vector9d pB{};
  Kokkos::RangePolicy<execution_space> policy_remote(std::size_t{0}, n_local);
  auto reducerB = Reduction::make_kokkos_reducer<Utils::Vector9d>(
      // NOLINTNEXTLINE(bugprone-exception-escape)
      [=](std::size_t const i, Utils::Vector9d &psum) noexcept {
        auto const gi = offset + i;
        auto const &pos_i = pm[gi].pos;
        auto const &m_i = pm[gi].m;
        for (auto j = offset + n_local; j < n_total; ++j) {
          auto const &pos_j = pm[j].pos;
          auto const &m_j = pm[j].m;
          auto const d0 = with_replicas ? (pos_i - pos_j)
                                        : box_geo.get_mi_vector(pos_i, pos_j);
          for (std::size_t s = 0; s < n_shifts; ++s) {
            auto const rn = d0 + shifts_ptr[s];
            psum += Utils::flatten(
                Utils::tensor_product(rn, pair_force(rn, m_i, m_j).f));
          }
        }
      },
      join_op);
  Kokkos::parallel_reduce("dds_pressure_remote", policy_local, reducerB, pB);

  return prefactor * (pA + pB);
}

/**
 * @brief Calculate total dipole field of each particle.
 *
 * This employs a parallel N-square loop over all particles.
 * Logically this is equivalent to the potential calculation
 * in @ref DipolarDirectSum::long_range_energy, which calculates
 * a naive N-square sum. The difference is summation range,
 * and the kernel calculates the dipole field rather than the energy.
 */
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void DipolarDirectSum::dipole_field_at_part_cpu() const {
  assert(not m_is_gpu);
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &box_l = box_geo.length();
  auto const particles = system.cell_structure->local_particles();
  /* collect particle data */
  auto [local_particles, all_posmom, reqs, offset_signed] =
      gather_particle_data(box_geo, particles);

  auto const ncut = get_n_cut(box_geo, n_replicas);
  auto const with_replicas = (ncut.norm2() > 0);
  auto const shifts = make_image_shifts(ncut, box_l);

  auto const offset = static_cast<std::size_t>(offset_signed);
  auto const n_local = local_particles.size();
  auto const n_total = all_posmom.size();

  auto const prefactor_local = prefactor;

  /* The field sweeps over all j, so every view slice is needed. Unlike the
   * force/energy kernels there is no local-only computation to overlap with,
   * so wait for the remote data first, then fill all slices. */
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  /* Raw pointer to the gathered AoS data; all slices are now populated, and it
   * outlives the fence. Safe to capture by value in the Kokkos [=] lambda. */
  auto const *pm = all_posmom.data();

  using execution_space = Kokkos::DefaultExecutionSpace;

  /* Raw pointers so the Kokkos lambdas do not capture std::vector by value. */
  auto *local_particles_ptr = local_particles.data();
  auto const *shifts_ptr = shifts.data();
  auto const n_shifts = shifts.size();

  Kokkos::RangePolicy<execution_space> policy(std::size_t{0}, n_local);
  Kokkos::parallel_for("dds_dipole_field", policy, [=](std::size_t const i) {
    auto const gi = offset + i;
    auto const &pos_i = pm[gi].pos;
    auto const &m_i = pm[gi].m;
    Utils::Vector3d u{};

    /* (a) self-image term over shifts[1..] (primary excluded) */
    for (std::size_t s = 1; s < n_shifts; ++s)
      u += dipole_field(shifts_ptr[s], m_i);

    /* Sweep over all j in [0, n_total), self-primary excluded by splitting
     * the range into [0, gi) and [gi + 1, n_total). */
    std::size_t const ranges[2][2] = {{std::size_t{0}, gi}, {gi + 1, n_total}};
    for (auto const &range : ranges) {
      auto const range_begin = range[0];
      auto const range_end = range[1];
      for (auto j = range_begin; j < range_end; ++j) {
        auto const &pos_j = pm[j].pos;
        auto const &m_j = pm[j].m;
        auto const d0 = with_replicas ? (pos_i - pos_j)
                                      : box_geo.get_mi_vector(pos_i, pos_j);
        for (std::size_t s = 0; s < n_shifts; ++s)
          u += dipole_field(d0 + shifts_ptr[s], m_j);
      }
    }
    local_particles_ptr[i]->dip_fld() = prefactor_local * u;
  });
  Kokkos::fence();
}
#endif // ESPRESSO_DIPOLE_FIELD_TRACKING

DipolarDirectSum::DipolarDirectSum(double prefactor, int n_replicas, bool gpu) {
  set_prefactor(prefactor);
  m_is_gpu = gpu;
  this->n_replicas = n_replicas;
  if (n_replicas < 0) {
    throw std::domain_error("Parameter 'n_replicas' must be >= 0");
  }
}

#endif // ESPRESSO_DIPOLES
