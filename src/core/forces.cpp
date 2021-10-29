
/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "EspressoSystemInterface.hpp"

#include "collision.hpp"
#include "comfixed_global.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/icc.hpp"
#include "electrostatics_magnetostatics/p3m_gpu.hpp"
#include "forcecap.hpp"
#include "forces_inline.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "immersed_boundaries.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "short_range_loop.hpp"
#include "virtual_sites.hpp"

#include <profiler/profiler.hpp>

#include <cassert>

ActorList forceActors;

void init_forces(const ParticleRange &particles, double time_step, double kT) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  /* The force initialization depends on the used thermostat and the
     thermodynamic ensemble */

#ifdef NPT
  npt_reset_instantaneous_virials();
#endif

  /* initialize forces with Langevin thermostat forces
     or zero depending on the thermostat
     set torque to zero for all and rescale quaternions
  */
  for (auto &p : particles) {
    p.f = init_local_particle_force(p, time_step, kT);
  }

  /* initialize ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  for (auto &p : cell_structure.ghost_particles()) {
    p.f = init_ghost_force(p);
  }
}

void init_forces_ghosts(const ParticleRange &particles) {
  for (auto &p : particles) {
    p.f = init_ghost_force(p);
  }
}

void force_calc(CellStructure &cell_structure, double time_step, double kT) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  espressoSystemInterface.update();

#ifdef COLLISION_DETECTION
  prepare_local_collision_queue();
#endif

  auto particles = cell_structure.local_particles();
  auto ghost_particles = cell_structure.ghost_particles();
#ifdef ELECTROSTATICS
  icc_iteration(particles, cell_structure.ghost_particles());
#endif
  init_forces(particles, time_step, kT);

  for (auto &forceActor : forceActors) {
    forceActor->computeForces(espressoSystemInterface);
#ifdef ROTATION
    forceActor->computeTorques(espressoSystemInterface);
#endif
  }

  calc_long_range_forces(particles);

#ifdef ELECTROSTATICS
  auto const coulomb_cutoff = Coulomb::cutoff(box_geo.length());
#else
  auto const coulomb_cutoff = INACTIVE_CUTOFF;
#endif

#ifdef DIPOLES
  auto const dipole_cutoff = Dipole::cutoff(box_geo.length());
#else
  auto const dipole_cutoff = INACTIVE_CUTOFF;
#endif

  short_range_loop(
      add_bonded_force,
      [](Particle &p1, Particle &p2, Distance const &d) {
        add_non_bonded_pair_force(p1, p2, d.vec21, sqrt(d.dist2), d.dist2);
#ifdef COLLISION_DETECTION
        if (collision_params.mode != COLLISION_MODE_OFF)
          detect_collision(p1, p2, d.dist2);
#endif
      },
      maximal_cutoff(), maximal_cutoff_bonded(),
      VerletCriterion{skin, interaction_range(), coulomb_cutoff, dipole_cutoff,
                      collision_detection_cutoff()});

  Constraints::constraints.add_forces(particles, get_sim_time());

  if (max_oif_objects) {
    // There are two global quantities that need to be evaluated:
    // object's surface and object's volume.
    for (int i = 0; i < max_oif_objects; i++) {
      auto const area_volume = calc_oif_global(i, cell_structure);
      if (fabs(area_volume[0]) < 1e-100 && fabs(area_volume[1]) < 1e-100)
        break;
      add_oif_global_forces(area_volume, i, cell_structure);
    }
  }

  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries.volume_conservation(cell_structure);

  lb_lbcoupling_calc_particle_lattice_ia(thermo_virtual, particles,
                                         ghost_particles, time_step);

#ifdef CUDA
  copy_forces_from_GPU(particles);
#endif

// VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
  virtual_sites()->back_transfer_forces_and_torques();
#endif

  // Communication Step: ghost forces
  cell_structure.ghosts_reduce_forces();

  // should be pretty late, since it needs to zero out the total force
  comfixed.apply(comm_cart, particles);

  // Needs to be the last one to be effective
  forcecap_cap(particles);

  // mark that forces are now up-to-date
  recalc_forces = false;
}

void calc_long_range_forces(const ParticleRange &particles) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::calc_long_range_force(particles);

#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* calculate k-space part of the magnetostatic interaction. */
  Dipole::calc_long_range_force(particles);
#endif /*ifdef DIPOLES */
}

#ifdef NPT
void npt_add_virial_force_contribution(const Utils::Vector3d &force,
                                       const Utils::Vector3d &d) {
  npt_add_virial_contribution(force, d);
}
#endif
