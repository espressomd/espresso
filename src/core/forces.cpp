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

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "constraints/Constraints.hpp"
#include "electrostatics/icc.hpp"
#include "electrostatics/p3m_gpu.hpp"
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
#include "short_range_loop.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/variant.hpp>

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#include <cassert>
#include <cmath>
#include <memory>
#include <span>
#include <variant>

/** External particle forces */
static ParticleForce external_force(Particle const &p) {
  ParticleForce f = {};

#ifdef EXTERNAL_FORCES
  f.f += p.ext_force();
#ifdef ROTATION
  f.torque += p.ext_torque();
#endif
#endif

#ifdef ENGINE
  // apply a swimming force in the direction of
  // the particle's orientation axis
  if (p.swimming().swimming and !p.swimming().is_engine_force_on_fluid) {
    f.f += p.swimming().f_swim * p.calc_director();
  }
#endif

  return f;
}

static void init_forces(ParticleRange const &particles,
                        ParticleRange const &ghost_particles) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  for (auto &p : particles) {
    p.force_and_torque() = external_force(p);
  }

  init_forces_ghosts(ghost_particles);
}

void init_forces_ghosts(ParticleRange const &particles) {
  for (auto &p : particles) {
    p.force_and_torque() = {};
  }
}

static void force_capping(ParticleRange const &particles, double force_cap) {
  if (force_cap > 0.) {
    auto const force_cap_sq = Utils::sqr(force_cap);
    for (auto &p : particles) {
      auto const force_sq = p.force().norm2();
      if (force_sq > force_cap_sq) {
        p.force() *= force_cap / std::sqrt(force_sq);
      }
    }
  }
}

void System::System::calculate_forces() {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
#ifdef CUDA
#ifdef CALIPER
  CALI_MARK_BEGIN("copy_particles_to_GPU");
#endif
  gpu.update();
#ifdef CALIPER
  CALI_MARK_END("copy_particles_to_GPU");
#endif
#endif // CUDA

#ifdef COLLISION_DETECTION
  collision_detection->clear_queue();
#endif
  bond_breakage->clear_queue();
  auto particles = cell_structure->local_particles();
  auto ghost_particles = cell_structure->ghost_particles();
#ifdef ELECTROSTATICS
  if (coulomb.impl->extension) {
    if (auto icc = std::get_if<std::shared_ptr<ICCStar>>(
            get_ptr(coulomb.impl->extension))) {
      (**icc).iteration(*cell_structure, particles, ghost_particles);
    }
  }
#endif // ELECTROSTATICS
#ifdef NPT
  npt_reset_instantaneous_virials();
#endif
  init_forces(particles, ghost_particles);
  thermostat_force_init();

  calc_long_range_forces(particles);

  auto const elc_kernel = coulomb.pair_force_elc_kernel();
  auto const coulomb_kernel = coulomb.pair_force_kernel();
  auto const dipoles_kernel = dipoles.pair_force_kernel();

#ifdef ELECTROSTATICS
  auto const coulomb_cutoff = coulomb.cutoff();
#else
  auto const coulomb_cutoff = INACTIVE_CUTOFF;
#endif

#ifdef DIPOLES
  auto const dipole_cutoff = dipoles.cutoff();
#else
  auto const dipole_cutoff = INACTIVE_CUTOFF;
#endif
#ifdef COLLISION_DETECTION
  auto const collision_detection_cutoff = collision_detection->cutoff();
#else
  auto const collision_detection_cutoff = INACTIVE_CUTOFF;
#endif

  short_range_loop(
      [coulomb_kernel_ptr = get_ptr(coulomb_kernel), &bonded_ias = *bonded_ias,
       &bond_breakage = *bond_breakage, &box_geo = *box_geo](
          Particle &p1, int bond_id, std::span<Particle *> partners) {
        return add_bonded_force(p1, bond_id, partners, bonded_ias,
                                bond_breakage, box_geo, coulomb_kernel_ptr);
      },
      [coulomb_kernel_ptr = get_ptr(coulomb_kernel),
       dipoles_kernel_ptr = get_ptr(dipoles_kernel),
       elc_kernel_ptr = get_ptr(elc_kernel), &nonbonded_ias = *nonbonded_ias,
       &thermostat = *thermostat, &bonded_ias = *bonded_ias,
#ifdef COLLISION_DETECTION
       &collision_detection = *collision_detection,
#endif
       &box_geo = *box_geo](Particle &p1, Particle &p2, Distance const &d) {
        auto const &ia_params =
            nonbonded_ias.get_ia_param(p1.type(), p2.type());
        add_non_bonded_pair_force(p1, p2, d.vec21, sqrt(d.dist2), d.dist2,
                                  ia_params, thermostat, box_geo, bonded_ias,
                                  coulomb_kernel_ptr, dipoles_kernel_ptr,
                                  elc_kernel_ptr);
#ifdef COLLISION_DETECTION
        if (collision_detection.mode != CollisionModeType::OFF) {
          collision_detection.detect_collision(p1, p2, d.dist2);
        }
#endif
      },
      *cell_structure, maximal_cutoff(), bonded_ias->maximal_cutoff(),
      VerletCriterion<>{*this, cell_structure->get_verlet_skin(),
                        get_interaction_range(), coulomb_cutoff, dipole_cutoff,
                        collision_detection_cutoff});

  constraints->add_forces(particles, get_sim_time());
  oif_global->calculate_forces();

  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries->volume_conservation(*cell_structure);

  if (thermostat->lb and (propagation->used_propagations &
                          PropagationMode::TRANS_LB_MOMENTUM_EXCHANGE)) {
    lb_couple_particles();
  }

#ifdef CUDA
#ifdef CALIPER
  CALI_MARK_BEGIN("copy_forces_from_GPU");
#endif
  gpu.copy_forces_to_host(particles, this_node);
#ifdef CALIPER
  CALI_MARK_END("copy_forces_from_GPU");
#endif
#endif // CUDA

#ifdef VIRTUAL_SITES_RELATIVE
  if (propagation->used_propagations &
      (PropagationMode::TRANS_VS_RELATIVE | PropagationMode::ROT_VS_RELATIVE)) {
    vs_relative_back_transfer_forces_and_torques(*cell_structure);
  }
#endif

  // Communication step: ghost forces
  cell_structure->ghosts_reduce_forces();

  // should be pretty late, since it needs to zero out the total force
  comfixed->apply(particles);

  // Needs to be the last one to be effective
  force_capping(particles, force_cap);

  // mark that forces are now up-to-date
  propagation->recalc_forces = false;
}

void calc_long_range_forces(const ParticleRange &particles) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::get_coulomb().calc_long_range_force(particles);
#endif // ELECTROSTATICS

#ifdef DIPOLES
  /* calculate k-space part of the magnetostatic interaction. */
  Dipoles::get_dipoles().calc_long_range_force(particles);
#endif // DIPOLES
}

#ifdef NPT
void npt_add_virial_force_contribution(const Utils::Vector3d &force,
                                       const Utils::Vector3d &d) {
  npt_add_virial_contribution(force, d);
}
#endif
