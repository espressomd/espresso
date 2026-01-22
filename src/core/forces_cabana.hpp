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

#include <config/config.hpp>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM

#include "aosoa_pack.hpp"
#include "forces_inline.hpp"

#include <utils/Vector.hpp>

#include <Cabana_Core.hpp>

#include <omp.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <variant>
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
#ifdef ESPRESSO_P3M
  CoulombP3M const *p3m;
#endif

  ForcesKernel(
      BondedInteractionsMap const &bonded_ias_,
      InteractionsNonBonded const &nonbonded_ias_,
      Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel_,
      Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel_,
      Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel_,
      Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_u_kernel_,
      Coulomb::Solver const &coulomb_,
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
#ifdef ESPRESSO_P3M
    p3m = nullptr;
    if (auto &solver = coulomb_.impl->solver; solver.has_value()) {
      if (std::holds_alternative<std::shared_ptr<CoulombP3M>>(*solver)) {
        p3m = std::get<std::shared_ptr<CoulombP3M>>(*solver).get();
      }
    }
#endif
  }

  // Helper functions to check if specific algorithms are active
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
  gay_berne_active(double dist, IA_parameters const &ia_params) const {
#ifdef ESPRESSO_GAY_BERNE
    return dist < ia_params.gay_berne.cut;
#else
    return false;
#endif
  }

#ifdef ESPRESSO_NPT
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool npt_active() const {
    return global_virial != nullptr;
  }
#endif

#ifdef ESPRESSO_THOLE
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
  thole_active(IA_parameters const &ia_params) const {
    return (ia_params.thole.scaling_coeff != 0. and
            ia_params.thole.q1q2 != 0. and coulomb_kernel != nullptr);
  }
#endif

#ifdef ESPRESSO_DIPOLES
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool
  dipoles_active() const {
    return dipoles_kernel != nullptr;
  }
#endif

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t i, std::size_t j) const {

    auto const &ia_params =
        nonbonded_ias.get_ia_param(aosoa.type(i), aosoa.type(j));

    ParticleForce pf{};

    // calc distance
    auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
    auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
    auto const d = box_geo.get_mi_vector(pos1, pos2);
    auto const dist = d.norm();

    // Determine which data needs to be loaded based on active algorithms
#if defined(ESPRESSO_DIPOLES) or defined(ESPRESSO_GAY_BERNE)
    bool const need_directors =
        gay_berne_active(dist, ia_params) or dipoles_active();
#endif
#if defined(ESPRESSO_EXCLUSIONS) or defined(ESPRESSO_THOLE)
    bool const need_particle_pointers = aosoa.has_exclusion(i) or
                                        aosoa.has_exclusion(j) or
                                        thole_active(ia_params);
    Particle const *p1_ptr = nullptr;
    Particle const *p2_ptr = nullptr;
    if (need_particle_pointers) {
      p1_ptr = unique_particles.at(i);
      p2_ptr = unique_particles.at(j);
    }
#endif

    // Load directors only if needed
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
    Utils::Vector3d dir1{}, dir2{};
    if (need_directors) {
      dir1 = aosoa.get_vector_at(aosoa.director, i);
      dir2 = aosoa.get_vector_at(aosoa.director, j);
    }
#endif

    /***********************************************/
    /* non-bonded pair potentials                  */
    /***********************************************/

    if (dist < ia_params.max_cut) {
#ifdef ESPRESSO_EXCLUSIONS
      bool skip_non_bonded = false;
      if (aosoa.has_exclusion(i) or aosoa.has_exclusion(j)) {
        skip_non_bonded = not do_nonbonded(*p1_ptr, *p2_ptr);
      }
#else
      constexpr bool skip_non_bonded = false;
#endif
      if (not skip_non_bonded) {
        pf.f += calc_central_radial_force(ia_params, d, dist);

        // Only call Thole force kernel if active
#ifdef ESPRESSO_THOLE
        if (thole_active(ia_params)) {
          pf.f += thole_pair_force(*p1_ptr, *p2_ptr, ia_params, d, dist,
                                   bonded_ias, coulomb_kernel);
        }
#endif
        // Only call Gay-Berne force kernel if active
#ifdef ESPRESSO_GAY_BERNE
        if (gay_berne_active(dist, ia_params)) {
          pf += gb_pair_force(dir1, dir2, ia_params, d, dist);
        }
#endif
      } // not skip_non_bonded
    }

    /*********************************************************************/
    /* everything before this contributes to the virial pressure in NpT, */
    /* but nothing afterwards, since the contribution to pressure from   */
    /* electrostatic is calculated by energy                             */
    /*********************************************************************/
#ifdef ESPRESSO_NPT
    Utils::Vector3d virial{};
    if (npt_active()) {
      virial = hadamard_product(pf.f, d);
    }
#endif // ESPRESSO_NPT

    /***********************************************/
    /* thermostat                                  */
    /***********************************************/

    /* The inter dpd force should not be part of the virial */
#ifdef ESPRESSO_DPD
    if (thermostat.thermo_switch & THERMO_DPD) {
      auto const dist2 = dist * dist;
      auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
      auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
      auto const force =
          dpd_pair_force(pos1, vel1, aosoa.id(i), pos2, vel2, aosoa.id(j),
                         *thermostat.dpd, box_geo, ia_params, d, dist, dist2);
      pf += force;
    }
#endif // ESPRESSO_DPD

#ifdef ESPRESSO_ELECTROSTATICS
    Utils::Vector3d f1_asym{};
    Utils::Vector3d f2_asym{};
    // real-space electrostatic charge-charge interaction
    if (coulomb_kernel != nullptr) {
      auto const q1q2 = aosoa.charge(i) * aosoa.charge(j);
      if (q1q2 != 0) {
#ifdef ESPRESSO_P3M
        if (p3m) [[likely]] {
          pf.f += p3m->pair_force(q1q2, d, dist);
        } else
#endif
        {
          pf.f += (*coulomb_kernel)(q1q2, d, dist);
        }
        if (elc_kernel) {
          (*elc_kernel)(pos1, pos2, f1_asym, f2_asym, q1q2);
        }
#ifdef ESPRESSO_NPT
        if (npt_active()) {
          virial[0] += (*coulomb_u_kernel)(pos1, pos2, q1q2, d, dist);
        }
#endif // ESPRESSO_NPT
      }
    }
#endif // ESPRESSO_ELECTROSTATICS

    // Only call dipole force kernel if active
#ifdef ESPRESSO_DIPOLES
    if (dipoles_active()) {
      auto const d1d2 = aosoa.dipm(i) * aosoa.dipm(j);
      if (d1d2 != 0.) {
        pf += (*dipoles_kernel)(d1d2, aosoa.dipm(i) * dir1,
                                aosoa.dipm(j) * dir2, d, dist, dist * dist);
      }
    }
#endif // ESPRESSO_DIPOLES

    auto opf = calc_opposing_force(pf, d);
#ifdef ESPRESSO_ELECTROSTATICS
    pf.f += f1_asym;
    opf.f += f2_asym;
#endif // ESPRESSO_ELECTROSTATICS

    auto const thread_id = omp_get_thread_num();

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
    if (npt_active()) {
      local_virial(thread_id, 0) += virial[0];
      local_virial(thread_id, 1) += virial[1];
      local_virial(thread_id, 2) += virial[2];
    }
#endif
  }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
