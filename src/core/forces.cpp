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
/** \file
 *  Force calculation.
 *
 *  The corresponding header file is forces.hpp.
 */

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "collision_detection/CollisionDetection.hpp"
#include "communication.hpp"
#include "constraints/Constraints.hpp"
#include "electrostatics/icc.hpp"
#include "forces_inline.hpp"
#include "galilei/ComFixed.hpp"
#include "immersed_boundary/ImmersedBoundaries.hpp"
#include "integrators/Propagation.hpp"
#include "lb/particle_coupling.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rotation.hpp"
#include "short_range_cabana.hpp"
#include "short_range_loop.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Cabana_Core.hpp>
#endif

#include <cassert>
#include <cmath>
#include <memory>
#include <span>
#include <variant>

/** External particle forces */
static ParticleForce external_force(Particle const &p) {
  ParticleForce f = {};

#ifdef ESPRESSO_EXTERNAL_FORCES
  f.f += p.ext_force();
#ifdef ESPRESSO_ROTATION
  f.torque += p.ext_torque();
#endif
#endif

#ifdef ESPRESSO_ENGINE
  // apply a swimming force in the direction of
  // the particle's orientation axis
  if (p.swimming().swimming and !p.swimming().is_engine_force_on_fluid) {
    f.f += p.swimming().f_swim * p.calc_director();
  }
#endif

  return f;
}

/** Combined force initialization and Langevin noise application */
void init_forces_and_thermostat(System::System const &system) {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  auto &cell_structure = *system.cell_structure;
  auto const &propagation = *system.propagation;
  auto const &thermostat = *system.thermostat;
  auto const kT = thermostat.kT;
  auto const time_step = system.get_time_step();

  // Check if Langevin thermostat is active
  bool const langevin_active =
      thermostat.langevin &&
      (propagation.used_propagations &
       (PropagationMode::TRANS_LANGEVIN | PropagationMode::ROT_LANGEVIN));

  // Single pass over all local particles
  cell_structure.for_each_local_particle([&](Particle &p) {
    // Initialize force with external forces
    p.force_and_torque() = external_force(p);

    // Apply Langevin noise if thermostat is active
    if (langevin_active) {
      auto const &langevin = *thermostat.langevin;
      if (propagation.should_propagate_with(p, PropagationMode::TRANS_LANGEVIN))
        p.force() += friction_thermo_langevin(langevin, p, time_step, kT);
#ifdef ESPRESSO_ROTATION
      if (propagation.should_propagate_with(p, PropagationMode::ROT_LANGEVIN))
        p.torque() += convert_vector_body_to_space(
            p, friction_thermo_langevin_rotation(langevin, p, time_step, kT));
#endif
    }
  });
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  cell_structure.reset_local_force();
#endif

  // Initialize ghost forces (unchanged)
  init_forces_ghosts(cell_structure);
}

void init_forces_ghosts(const CellStructure &cell_structure) {
  cell_structure.for_each_ghost_particle(
      [](Particle &p) { p.force_and_torque() = {}; });
}

static void force_capping(CellStructure &cell_structure, double force_cap) {
  if (force_cap > 0.) {
    auto const force_cap_sq = Utils::sqr(force_cap);
    cell_structure.for_each_local_particle(
        [&force_cap, &force_cap_sq](Particle &p) {
          auto const force_sq = p.force().norm2();
          if (force_sq > force_cap_sq) {
            p.force() *= force_cap / std::sqrt(force_sq);
          }
        });
  }
}

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
static void reinit_dip_fld(CellStructure const &cell_structure) {
  cell_structure.for_each_local_particle(
      [](Particle &p) { p.dip_fld() = {0., 0., 0.}; });
}
#endif

void System::System::calculate_forces() {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
#ifdef ESPRESSO_CUDA
  {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("copy_particles_to_GPU");
#endif
    gpu.update();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("copy_particles_to_GPU");
#endif
  }
#endif // ESPRESSO_CUDA

#ifdef ESPRESSO_COLLISION_DETECTION
  collision_detection->clear_queue();
  auto const collision_detection_cutoff = collision_detection->cutoff();
#else
  auto const collision_detection_cutoff = inactive_cutoff;
#endif
  bond_breakage->clear_queue();
  auto particles = cell_structure->local_particles();
#ifdef ESPRESSO_NPT
  if (propagation->used_propagations & PropagationMode::TRANS_LANGEVIN_NPT) {
    // reset virial part of instantaneous pressure
    npt_inst_pressure->p_vir = Utils::Vector3d{};
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  // reset dipole field
  reinit_dip_fld(*cell_structure);
#endif

  // Use combined function instead of two separate calls

  auto const elc_kernel = coulomb.pair_force_elc_kernel();
  auto const coulomb_kernel = coulomb.pair_force_kernel();
  auto const dipoles_kernel = dipoles.pair_force_kernel();
  auto const coulomb_u_kernel = coulomb.pair_energy_kernel();
  auto *const virial = get_npt_virial();

  // interaction kernel is defined
  auto bond_kernel = [coulomb_kernel_ptr = get_ptr(coulomb_kernel),
                      &bonded_ias = *bonded_ias,
                      &bond_breakage = *bond_breakage, virial,
                      &box_geo = *box_geo](Particle &p1, int bond_id,
                                           std::span<Particle *> partners) {
    return add_bonded_force(p1, bond_id, partners, bonded_ias, bond_breakage,
                            box_geo, virial, coulomb_kernel_ptr);
  };

  VerletCriterion<> const verlet_criterion{*this,
                                           cell_structure->get_verlet_skin(),
                                           get_interaction_range(),
                                           coulomb.cutoff(),
                                           dipoles.cutoff(),
                                           collision_detection_cutoff};

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  update_cabana_state(*cell_structure, verlet_criterion,
                      get_interaction_range(), propagation->integ_switch);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  if (coulomb.impl->extension) {
    update_icc_particles();
  }
#endif // ESPRESSO_ELECTROSTATICS
  init_forces_and_thermostat(*this);
  calc_long_range_forces(particles);

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#ifdef ESPRESSO_CALIPER
  CALI_MARK_BEGIN("cabana_short_range");
#endif
  using execution_space = Kokkos::DefaultExecutionSpace;
  auto const &unique_particles = cell_structure->get_unique_particles();
  auto const &local_force = cell_structure->get_local_force();
#ifdef ESPRESSO_ROTATION
  auto const &local_torque = cell_structure->get_local_torque();
#endif
#ifdef ESPRESSO_NPT
  auto const &local_virial = cell_structure->get_local_virial();
#endif
  auto const &aosoa = cell_structure->get_aosoa();

  ForcesKernel first_neighbor_kernel(
      *bonded_ias, *nonbonded_ias, get_ptr(coulomb_kernel),
      get_ptr(dipoles_kernel), get_ptr(elc_kernel), get_ptr(coulomb_u_kernel),
      *thermostat, *box_geo, unique_particles, local_force,
#ifdef ESPRESSO_ROTATION
      local_torque,
#endif
#ifdef ESPRESSO_NPT
      virial, local_virial,
#endif
      aosoa);

  cabana_short_range(bond_kernel, first_neighbor_kernel, *cell_structure,
                     get_interaction_range(), bonded_ias->maximal_cutoff(),
                     verlet_criterion, propagation->integ_switch);
  // Force and Torque reduction
  int num_threads = execution_space().concurrency();
  Kokkos::RangePolicy<execution_space> policy(std::size_t{0},
                                              unique_particles.size());
  Kokkos::parallel_for("reduction", policy,
                       [&local_force,
#ifdef ESPRESSO_ROTATION
                        &local_torque,
#endif
                        &unique_particles, num_threads](std::size_t const i) {
                         Utils::Vector3d force{};
#ifdef ESPRESSO_ROTATION
                         Utils::Vector3d torque{};
#endif
                         for (int tid = 0; tid < num_threads; ++tid) {
                           force[0] += local_force(i, tid, 0);
                           force[1] += local_force(i, tid, 1);
                           force[2] += local_force(i, tid, 2);
#ifdef ESPRESSO_ROTATION
                           torque[0] += local_torque(i, tid, 0);
                           torque[1] += local_torque(i, tid, 1);
                           torque[2] += local_torque(i, tid, 2);
#endif
                         }
                         unique_particles.at(i)->force() += force;
#ifdef ESPRESSO_ROTATION
                         unique_particles.at(i)->torque() += torque;
#endif
                       });
  Kokkos::fence();

#ifdef ESPRESSO_NPT
  if (virial) {
    for (int tid = 0; tid < num_threads; ++tid) {
      (*virial)[0] += local_virial(tid, 0);
      (*virial)[1] += local_virial(tid, 1);
      (*virial)[2] += local_virial(tid, 2);
    }
  }
#endif

#ifdef ESPRESSO_COLLISION_DETECTION
  auto collision_kernel = [&collision_detection = *collision_detection](
                              Particle const &p1, Particle const &p2,
                              Distance const &d) {
    collision_detection.detect_collision(p1, p2, d.dist2);
  };
  if (not collision_detection->is_off()) {
    cell_structure->non_bonded_loop(collision_kernel, verlet_criterion);
  }
#endif

#ifdef ESPRESSO_CALIPER
  CALI_MARK_END("cabana_short_range");
#endif

#else // ESPRESSO_SHARED_MEMORY_PARALLELISM

#ifdef ESPRESSO_CALIPER
  CALI_MARK_BEGIN("serial_short_range");
#endif

  auto pair_kernel = [coulomb_kernel_ptr = get_ptr(coulomb_kernel),
                      dipoles_kernel_ptr = get_ptr(dipoles_kernel),
                      elc_kernel_ptr = get_ptr(elc_kernel),
                      coulomb_u_kernel_ptr = get_ptr(coulomb_u_kernel),
                      &nonbonded_ias = *nonbonded_ias,
                      &thermostat = *thermostat, &bonded_ias = *bonded_ias,
                      virial,
#ifdef ESPRESSO_COLLISION_DETECTION
                      &collision_detection = *collision_detection,
#endif
                      &box_geo = *box_geo](Particle &p1, Particle &p2,
                                           Distance const &d) {
    auto const &ia_params = nonbonded_ias.get_ia_param(p1.type(), p2.type());
    add_non_bonded_pair_force(
        p1, p2, d.vec21, sqrt(d.dist2), d.dist2, p1.q() * p2.q(), ia_params,
        thermostat, box_geo, bonded_ias, virial, coulomb_kernel_ptr,
        dipoles_kernel_ptr, elc_kernel_ptr, coulomb_u_kernel_ptr);
#ifdef ESPRESSO_COLLISION_DETECTION
    if (not collision_detection.is_off()) {
      collision_detection.detect_collision(p1, p2, d.dist2);
    }
#endif
  };

  short_range_loop(bond_kernel, pair_kernel, *cell_structure, maximal_cutoff(),
                   bonded_ias->maximal_cutoff(), verlet_criterion);

#ifdef ESPRESSO_CALIPER
  CALI_MARK_END("serial_short_range");
#endif

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

  constraints->add_forces(particles, get_sim_time());
  oif_global->calculate_forces();

  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries->volume_conservation(*cell_structure);

  if (thermostat->lb and (propagation->used_propagations &
                          PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE)) {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("lb_particle_coupling");
#endif
    lb_couple_particles();
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("lb_particle_coupling");
#endif
  }

#ifdef ESPRESSO_CUDA
  {
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN("copy_forces_from_GPU");
#endif
    gpu.copy_forces_to_host(particles, this_node);

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    gpu.copy_dip_fld_to_host(particles, this_node);
#endif

#ifdef ESPRESSO_CALIPER
    CALI_MARK_END("copy_forces_from_GPU");
#endif
  }
#endif // ESPRESSO_CUDA

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  if (propagation->used_propagations &
      (PropagationMode::TRANS_VS_RELATIVE | PropagationMode::ROT_VS_RELATIVE |
       PropagationMode::ROT_VS_INDEPENDENT)) {
    vs_relative_back_transfer_forces_and_torques(*cell_structure);
  }
#endif

  // Communication step: ghost forces
  cell_structure->ghosts_reduce_forces();

  // should be pretty late, since it needs to zero out the total force
  comfixed->apply(particles);

  // Needs to be the last one to be effective
  force_capping(*cell_structure, force_cap);

  // mark that forces are now up-to-date
  propagation->recalc_forces = false;
}

void calc_long_range_forces(const ParticleRange &particles) {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

#ifdef ESPRESSO_ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::get_coulomb().calc_long_range_force(particles);
#endif // ESPRESSO_ELECTROSTATICS

#ifdef ESPRESSO_DIPOLES
  /* calculate k-space part of the magnetostatic interaction. */
  Dipoles::get_dipoles().calc_long_range_force(particles);
#endif // ESPRESSO_DIPOLES
}
