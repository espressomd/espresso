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

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM

#include "aosoa_pack.hpp"
#include "forces_inline.hpp"

#include <utils/Vector.hpp>

#include <Cabana_Core.hpp>

#include <omp.h>

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
  Dipoles::ShortRangeForceKernel::kernel_type const *const dipoles_kernel;
  Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel;
  Coulomb::ShortRangeEnergyKernel::kernel_type const *const coulomb_u_kernel;
  Thermostat::Thermostat const &thermostat;
  BoxGeometry const &box_geo;
  std::vector<Particle *> const &unique_particles;
  CellStructure::ForceType const &local_force;
#ifdef ESPRESSO_ROTATION
  CellStructure::ForceType const &local_torque;
#endif
#ifdef ESPRESSO_NPT
  Utils::Vector3d *const global_virial;
  CellStructure::VirialType const &local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;

  ForcesKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
      Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel_,
      Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
      Thermostat::Thermostat const &thermostat_, BoxGeometry const &box_geo_,
      std::vector<Particle *> const &unique_particles_,
      CellStructure::ForceType const &local_force_,
#ifdef ESPRESSO_ROTATION
      CellStructure::ForceType const &local_torque_,
#endif
#ifdef ESPRESSO_NPT
      Utils::Vector3d *const global_virial_,
      CellStructure::VirialType const &local_virial_,
#endif
      CellStructure::AoSoA_pack const &aosoa_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb_kernel(coulomb_kernel_), dipoles_kernel(dipoles_kernel_),
        elc_kernel(elc_kernel_), coulomb_u_kernel(coulomb_u_kernel_),
        thermostat(thermostat_), box_geo(box_geo_),
        unique_particles(unique_particles_), local_force(local_force_),
#ifdef ESPRESSO_ROTATION
        local_torque(local_torque_),
#endif
#ifdef ESPRESSO_NPT
        global_virial(global_virial_), local_virial(local_virial_),
#endif
        aosoa(aosoa_) {
  }

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t i, std::size_t j) const {

    auto const thread_id = omp_get_thread_num();

    auto const &ia_params =
        nonbonded_ias.get_ia_param(aosoa.type(i), aosoa.type(j));

    ParticleForce pf{};
    Utils::Vector3d const pos1 = {aosoa.position(i, 0), aosoa.position(i, 1),
                                  aosoa.position(i, 2)};
    Utils::Vector3d const pos2 = {aosoa.position(j, 0), aosoa.position(j, 1),
                                  aosoa.position(j, 2)};

#ifdef ESPRESSO_NPT
    Utils::Vector3d virial{};
    auto *const virial_handle = global_virial ? &virial : nullptr;
#endif
    auto const d = box_geo.get_mi_vector(pos1, pos2);
    auto const dist = d.norm();

#if defined(ESPRESSO_EXCLUSIONS) or defined(ESPRESSO_DPD) or                   \
    defined(ESPRESSO_DIPOLES)
    auto const &p1 = *unique_particles.at(i);
    auto const &p2 = *unique_particles.at(j);
#endif

    /***********************************************/
    /* non-bonded pair potentials                  */
    /***********************************************/

    if (dist < ia_params.max_cut) {
#ifdef ESPRESSO_EXCLUSIONS
      if (do_nonbonded(p1, p2)) {
#endif
        pf += calc_central_radial_force(ia_params, d, dist);
#ifdef ESPRESSO_THOLE
        pf.f += thole_pair_force(p1, p2, ia_params, d, dist, bonded_ias,
                                 coulomb_kernel);
#endif
#ifdef ESPRESSO_GAY_BERNE
        pf += calc_non_central_force(p1, p2, ia_params, d, dist);
#endif
#ifdef ESPRESSO_EXCLUSIONS
      }
#endif
    }

    /*********************************************************************/
    /* everything before this contributes to the virial pressure in NpT, */
    /* but nothing afterwards, since the contribution to pressure from   */
    /* electrostatic is calculated by energy                             */
    /*********************************************************************/
#ifdef ESPRESSO_NPT
    if (virial_handle) {
      *virial_handle += hadamard_product(pf.f, d);
    }
#endif // ESPRESSO_NPT

    /***********************************************/
    /* thermostat                                  */
    /***********************************************/

    /* The inter dpd force should not be part of the virial */
#ifdef ESPRESSO_DPD
    if (thermostat.thermo_switch & THERMO_DPD) {
      auto const dist2 = dist * dist;
      auto const force =
          dpd_pair_force(pos1, p1.v(), aosoa.id(i), pos2, p2.v(), aosoa.id(j),
                         *thermostat.dpd, box_geo, ia_params, d, dist, dist2);
      pf += force;
    }
#endif // ESPRESSO_DPD

#ifdef ESPRESSO_ELECTROSTATICS
    auto const q1q2 = aosoa.charge(i) * aosoa.charge(j);
    ParticleForce p1f_asym{};
    ParticleForce p2f_asym{};
    // real-space electrostatic charge-charge interaction
    if (q1q2 != 0. and coulomb_kernel != nullptr) {
      pf.f += (*coulomb_kernel)(q1q2, d, dist);
      if (elc_kernel) {
        (*elc_kernel)(pos1, pos2, p1f_asym, p2f_asym, q1q2);
      }
#ifdef ESPRESSO_NPT
      if (virial_handle) {
        (*virial_handle)[0] += (*coulomb_u_kernel)(pos1, pos2, q1q2, d, dist);
      }
#endif // ESPRESSO_NPT
    }
#endif // ESPRESSO_ELECTROSTATICS
#ifdef ESPRESSO_DIPOLES
    // real-space magnetic dipole-dipole interaction
    if (dipoles_kernel) {
      auto const d1d2 = p1.dipm() * p2.dipm();
      if (d1d2 != 0.) {
        pf += (*dipoles_kernel)(d1d2, p1.calc_dip(), p2.calc_dip(), d, dist,
                                dist * dist);
      }
    }
#endif // ESPRESSO_DIPOLES

    auto opf = calc_opposing_force(pf, d);
#ifdef ESPRESSO_ELECTROSTATICS
    pf += p1f_asym;
    opf += p2f_asym;
#endif // ESPRESSO_ELECTROSTATICS

    local_force(i, thread_id, 0) += pf.f[0];
    local_force(i, thread_id, 1) += pf.f[1];
    local_force(i, thread_id, 2) += pf.f[2];
#ifdef ESPRESSO_ROTATION
    local_torque(i, thread_id, 0) += pf.torque[0];
    local_torque(i, thread_id, 1) += pf.torque[1];
    local_torque(i, thread_id, 2) += pf.torque[2];
#endif

    local_force(j, thread_id, 0) += opf.f[0];
    local_force(j, thread_id, 1) += opf.f[1];
    local_force(j, thread_id, 2) += opf.f[2];
#ifdef ESPRESSO_ROTATION
    local_torque(j, thread_id, 0) += opf.torque[0];
    local_torque(j, thread_id, 1) += opf.torque[1];
    local_torque(j, thread_id, 2) += opf.torque[2];
#endif
#ifdef ESPRESSO_NPT
    if (virial_handle) {
      local_virial(thread_id, 0) += virial[0];
      local_virial(thread_id, 1) += virial[1];
      local_virial(thread_id, 2) += virial[2];
    }
#endif
  }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
