/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "aosoa_pack.hpp"
#include "forces_inline.hpp"
#include "short_range_cabana_helpers.hpp"

#include <utils/Vector.hpp>

#include <Cabana_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <cstddef>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

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
  CellStructure::ScatterForce local_force;
#ifdef ESPRESSO_ROTATION
  CellStructure::ScatterForce local_torque;
#endif
#ifdef ESPRESSO_NPT
  Utils::Vector3d *const global_virial;
  CellStructure::ScatterVirial local_virial;
#endif
  CellStructure::AoSoA_pack const &aosoa;
#ifdef ESPRESSO_P3M
  CoulombP3M const *p3m;
#endif
  double system_max_cutoff_sq;

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
      CellStructure::ScatterForce local_force_,
#ifdef ESPRESSO_ROTATION
      CellStructure::ScatterForce local_torque_,
#endif
#ifdef ESPRESSO_NPT
      Utils::Vector3d *const global_virial_,
      CellStructure::ScatterVirial local_virial_,
#endif
      CellStructure::AoSoA_pack const &aosoa_, double system_max_cutoff_)
      : bonded_ias(bonded_ias_), nonbonded_ias(nonbonded_ias_),
        coulomb_kernel(coulomb_kernel_), dipoles_kernel(dipoles_kernel_),
        elc_kernel(elc_kernel_), coulomb_u_kernel(coulomb_u_kernel_),
        thermostat(thermostat_), box_geo(box_geo_),
        unique_particles(unique_particles_),
        local_force(std::move(local_force_)),
#ifdef ESPRESSO_ROTATION
        local_torque(std::move(local_torque_)),
#endif
#ifdef ESPRESSO_NPT
        global_virial(global_virial_), local_virial(std::move(local_virial_)),
#endif
        aosoa(aosoa_), system_max_cutoff_sq(Utils::sqr(system_max_cutoff_)) {
#ifdef ESPRESSO_P3M
    p3m = nullptr;
    if (auto &solver = coulomb_.impl->solver; solver.has_value()) {
      if (std::holds_alternative<std::shared_ptr<CoulombP3M>>(*solver)) {
        p3m = std::get<std::shared_ptr<CoulombP3M>>(*solver).get();
      }
    }
#endif
  }

#ifdef ESPRESSO_NPT
  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION bool npt_active() const {
    return global_virial != nullptr;
  }
#endif

  ESPRESSO_ATTR_ALWAYS_INLINE KOKKOS_INLINE_FUNCTION void
  operator()(std::size_t i, std::size_t j) const {

    // calc distance (component-wise, avoids constructing pos1/pos2 Vector3d
    // on the hot early-exit path; pos1/pos2 are built lazily below only
    // where kernels actually require them)
    auto const d = box_geo.get_mi_vector(
        aosoa.position(i, 0), aosoa.position(i, 1), aosoa.position(i, 2),
        aosoa.position(j, 0), aosoa.position(j, 1), aosoa.position(j, 2));
    auto const dist_sq = d.norm2();

    // Early exit if distance > maximal global cutoff
    if (dist_sq > system_max_cutoff_sq)
      return;
    auto const dist = std::sqrt(dist_sq);
    auto const &ia_params =
        nonbonded_ias.get_ia_param(aosoa.type(i), aosoa.type(j));

    ParticleForce pf{};

    // Determine which data needs to be loaded based on active algorithms
#if defined(ESPRESSO_EXCLUSIONS) or defined(ESPRESSO_THOLE)
    bool need_particle_pointers = false;
#ifdef ESPRESSO_EXCLUSIONS
    need_particle_pointers |= aosoa.has_exclusion(i) or aosoa.has_exclusion(j);
#endif
#ifdef ESPRESSO_THOLE
    need_particle_pointers |=
        thole_active(ia_params, coulomb_kernel != nullptr);
#endif

    Particle const *p1_ptr = nullptr;
    Particle const *p2_ptr = nullptr;
    if (need_particle_pointers) {
      p1_ptr = unique_particles.at(i);
      p2_ptr = unique_particles.at(j);
    }
#endif

    /***********************************************/
    /* non-bonded pair potentials                  */
    /***********************************************/

    if (dist <= ia_params.max_cut) {
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
        if (thole_active(ia_params, coulomb_kernel != nullptr)) {
          pf.f += thole_pair_force(*p1_ptr, *p2_ptr, ia_params, d, dist,
                                   bonded_ias, coulomb_kernel);
        }
#endif
        // Only call Gay-Berne force kernel if active
#ifdef ESPRESSO_GAY_BERNE
        if (gay_berne_active(dist, ia_params)) {
          auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
          auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
          pf += gb_pair_force(dir1, dir2, ia_params, d, dist);
        }
#endif
      } // not skip_non_bonded
    } // not dist > ia_params.max_cut

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
    if (dpd_active(ia_params, thermostat.thermo_switch)) {
      auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
      auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
      auto const vel1 = aosoa.get_vector_at(aosoa.velocity, i);
      auto const vel2 = aosoa.get_vector_at(aosoa.velocity, j);
      auto const force =
          dpd_pair_force(pos1, vel1, aosoa.id(i), pos2, vel2, aosoa.id(j),
                         *thermostat.dpd, box_geo, ia_params, d, dist, dist_sq);
      pf += force;
    }
#endif // ESPRESSO_DPD

#ifdef ESPRESSO_ELECTROSTATICS
    Utils::Vector3d f1_asym{};
    Utils::Vector3d f2_asym{};
    // real-space electrostatic charge-charge interaction
    if (coulomb_kernel != nullptr) {
      if ((aosoa.charge(i) != 0.) and (aosoa.charge(j) != 0.)) {
        auto const q1q2 = aosoa.charge(i) * aosoa.charge(j);
#ifdef ESPRESSO_P3M
        if (p3m) [[likely]] {
          pf.f += p3m->pair_force(q1q2, d, dist);
        } else
#endif
        {
          pf.f += (*coulomb_kernel)(q1q2, d, dist);
        }
        if (elc_kernel) {
          auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
          auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
          (*elc_kernel)(pos1, pos2, f1_asym, f2_asym, q1q2);
        }
#ifdef ESPRESSO_NPT
        if (npt_active()) {
          auto const pos1 = aosoa.get_vector_at(aosoa.position, i);
          auto const pos2 = aosoa.get_vector_at(aosoa.position, j);
          virial[0] += (*coulomb_u_kernel)(pos1, pos2, q1q2, d, dist);
        }
#endif // ESPRESSO_NPT
      }
    }
#endif // ESPRESSO_ELECTROSTATICS

    // Only call dipole force kernel if active
#ifdef ESPRESSO_DIPOLES
    if (dipoles_kernel != nullptr) {
      auto const d1d2 = aosoa.dipm(i) * aosoa.dipm(j);
      if (d1d2 != 0.) {
        auto const dir1 = aosoa.get_vector_at(aosoa.director, i);
        auto const dir2 = aosoa.get_vector_at(aosoa.director, j);
        pf += (*dipoles_kernel)(d1d2, aosoa.dipm(i) * dir1,
                                aosoa.dipm(j) * dir2, d, dist, dist_sq);
      }
    }
#endif // ESPRESSO_DIPOLES

    auto opf = calc_opposing_force(pf, d);
#ifdef ESPRESSO_ELECTROSTATICS
    pf.f += f1_asym;
    opf.f += f2_asym;
#endif // ESPRESSO_ELECTROSTATICS

    auto access_force = local_force.access();

    access_force(i, 0) += pf.f[0];
    access_force(i, 1) += pf.f[1];
    access_force(i, 2) += pf.f[2];
#ifdef ESPRESSO_ROTATION
    auto access_torque = local_torque.access();
    access_torque(i, 0) += pf.torque[0];
    access_torque(i, 1) += pf.torque[1];
    access_torque(i, 2) += pf.torque[2];
#endif

    access_force(j, 0) += opf.f[0];
    access_force(j, 1) += opf.f[1];
    access_force(j, 2) += opf.f[2];
#ifdef ESPRESSO_ROTATION
    access_torque(j, 0) += opf.torque[0];
    access_torque(j, 1) += opf.torque[1];
    access_torque(j, 2) += opf.torque[2];
#endif
#ifdef ESPRESSO_NPT
    if (npt_active()) {
      auto access_virial = local_virial.access();
      access_virial(0) += virial[0];
      access_virial(1) += virial[1];
      access_virial(2) += virial[2];
    }
#endif
  }
};
