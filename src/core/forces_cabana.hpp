/*
 * Copyright (C) 2025 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#ifdef SHARED_MEMORY_PARALLELISM

#include "aosoa_pack.hpp"
#include "forces_inline.hpp"

#include <Cabana_Core.hpp>

#include <vector>

#if defined(__GNUG__) or defined(__clang__)
#define ESPRESSO_ATTR_ALWAYS_INLINE [[gnu::always_inline]]
#else
#define ESPRESSO_ATTR_ALWAYS_INLINE
#endif

struct ForcesKernel {
  BondedInteractionsMap const &bonded_ias;
  InteractionsNonBonded const &nonbonded_ias;
  Coulomb::ShortRangeForceKernel::kernel_type const *const coulomb_kernel;
#if defined(LONG_RANGE_KERNELS)
  Dipoles::ShortRangeForceKernel::kernel_type const *const dipoles_kernel;
  Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel;
  Coulomb::ShortRangeEnergyKernel::kernel_type const *const coulomb_u_kernel;
  Thermostat::Thermostat const &thermostat;
#endif
  BoxGeometry const &box_geo;
#if defined(LONG_RANGE_KERNELS) or defined(EXCLUSIONS)
  std::vector<Particle *> const &unique_particles;
#endif
  CellStructure::ForceType const &local_force;
#ifdef ROTATION
  CellStructure::ForceType const &local_torque;
#endif
#ifdef NPT
  CellStructure::VirialType const &local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;

  ForcesKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
#if defined(LONG_RANGE_KERNELS)
      Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel_,
      Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
      Thermostat::Thermostat const &thermostat_,
#endif
      BoxGeometry const &box_geo_,
#if defined(LONG_RANGE_KERNELS) or defined(EXCLUSIONS)
      std::vector<Particle *> const &unique_particles_,
#endif
      CellStructure::ForceType const &local_force_,
#ifdef ROTATION
      CellStructure::ForceType const &local_torque_,
#endif
#ifdef NPT
      CellStructure::VirialType const &local_virial_,
#endif
      CellStructure::AoSoA_pack const &aosoa_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb_kernel(coulomb_kernel_),
#if defined(LONG_RANGE_KERNELS)
        dipoles_kernel(dipoles_kernel_), elc_kernel(elc_kernel_),
        coulomb_u_kernel(coulomb_u_kernel_), thermostat(thermostat_),
#endif
        box_geo(box_geo_),
#if defined(LONG_RANGE_KERNELS) or defined(EXCLUSIONS)
        unique_particles(unique_particles_),
#endif
        local_force(local_force_),
#ifdef ROTATION
        local_torque(local_torque_),
#endif
#ifdef NPT
        local_virial(local_virial_),
#endif
        aosoa(aosoa_) {
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(int i, int j) const {

    auto const thread_id = omp_get_thread_num();

    auto const &ia_params =
        nonbonded_ias.get_ia_param(aosoa.type(i), aosoa.type(j));

    ParticleForce pf{};
#ifdef NPT
    Utils::Vector3d virial{};
#endif
    auto const d = box_geo.get_mi_vector(
        aosoa.position(i, 0), aosoa.position(i, 1), aosoa.position(i, 2),
        aosoa.position(j, 0), aosoa.position(j, 1), aosoa.position(j, 2));
    auto const dist = d.norm();

#if defined(LONG_RANGE_KERNELS) or defined(EXCLUSIONS)
    auto &p1 = *unique_particles.at(i);
    auto &p2 = *unique_particles.at(j);
#endif

#ifdef EXCLUSIONS
    auto const do_nonbonded_flag = do_nonbonded(p1, p2);
#else
#if defined(LONG_RANGE_KERNELS)
    auto constexpr do_nonbonded_flag = true;
#endif // LONG_RANGE_KERNELS
#endif

    if (dist < ia_params.max_cut) {
#ifdef EXCLUSIONS
      if (do_nonbonded_flag) {
#endif
        pf += calc_central_radial_force(ia_params, d, dist);
#ifdef EXCLUSIONS
      }
#endif
    }

#if defined(LONG_RANGE_KERNELS)
#ifdef ELECTROSTATICS
    auto const q1q2 = aosoa.charge(i) * aosoa.charge(j);
#else
    auto constexpr q1q2 = 0.;
#endif
    add_non_bonded_pair_force_with_p(p1, p2, pf,
#ifdef NPT
                                     virial,
#endif // NPT
                                     d, dist, dist * dist, q1q2, ia_params,
                                     do_nonbonded_flag, thermostat, box_geo,
                                     bonded_ias, coulomb_kernel, dipoles_kernel,
                                     elc_kernel, coulomb_u_kernel);
#endif // LONG_RANGE_KERNELS

    local_force(i, thread_id, 0) += pf.f[0];
    local_force(i, thread_id, 1) += pf.f[1];
    local_force(i, thread_id, 2) += pf.f[2];
#ifdef ROTATION
    local_torque(i, thread_id, 0) += pf.torque[0];
    local_torque(i, thread_id, 1) += pf.torque[1];
    local_torque(i, thread_id, 2) += pf.torque[2];
#endif

    auto const opf = calc_opposing_force(pf, d);
    local_force(j, thread_id, 0) += opf.f[0];
    local_force(j, thread_id, 1) += opf.f[1];
    local_force(j, thread_id, 2) += opf.f[2];
#ifdef ROTATION
    local_torque(j, thread_id, 0) += opf.torque[0];
    local_torque(j, thread_id, 1) += opf.torque[1];
    local_torque(j, thread_id, 2) += opf.torque[2];
#endif
#ifdef NPT
    local_virial(thread_id, 0) += virial[0];
    local_virial(thread_id, 1) += virial[1];
    local_virial(thread_id, 2) += virial[2];
#endif
  }
};

#endif // SHARED_MEMORY_PARALLELISM
