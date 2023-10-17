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
#include "bond_breakage/bond_breakage.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "electrostatics/icc.hpp"
#include "electrostatics/p3m_gpu.hpp"
#include "forcecap.hpp"
#include "forces_inline.hpp"
#include "galilei/ComFixed.hpp"
#include "immersed_boundaries.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "lb/particle_coupling.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "rotation.hpp"
#include "short_range_loop.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "thermostats/langevin_inline.hpp"
#include "virtual_sites.hpp"

#include <boost/variant.hpp>

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#include <cassert>
#include <memory>
#include <variant>

std::shared_ptr<ComFixed> comfixed = std::make_shared<ComFixed>();

/** Initialize the forces for a ghost particle */
inline ParticleForce init_ghost_force(Particle const &) { return {}; }

/** External particle forces */
inline ParticleForce external_force(Particle const &p) {
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

inline ParticleForce thermostat_force(Particle const &p, double time_step,
                                      double kT) {
  extern LangevinThermostat langevin;
  if (!(thermo_switch & THERMO_LANGEVIN)) {
    return {};
  }

#ifdef ROTATION
  return {friction_thermo_langevin(langevin, p, time_step, kT),
          p.can_rotate() ? convert_vector_body_to_space(
                               p, friction_thermo_langevin_rotation(
                                      langevin, p, time_step, kT))
                         : Utils::Vector3d{}};
#else
  return friction_thermo_langevin(langevin, p, time_step, kT);
#endif
}

/** Initialize the forces for a real particle */
inline ParticleForce init_real_particle_force(Particle const &p,
                                              double time_step, double kT) {
  return thermostat_force(p, time_step, kT) + external_force(p);
}

static void init_forces(const ParticleRange &particles,
                        const ParticleRange &ghost_particles, double time_step,
                        double kT) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

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
    p.force_and_torque() = init_real_particle_force(p, time_step, kT);
  }

  /* initialize ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  for (auto &p : ghost_particles) {
    p.force_and_torque() = init_ghost_force(p);
  }
}

void init_forces_ghosts(const ParticleRange &particles) {
  for (auto &p : particles) {
    p.force_and_torque() = init_ghost_force(p);
  }
}

void force_calc(CellStructure &cell_structure, double time_step, double kT) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  auto &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto const &nonbonded_ias = *system.nonbonded_ias;
#ifdef CUDA
#ifdef CALIPER
  CALI_MARK_BEGIN("copy_particles_to_GPU");
#endif
  system.gpu.update();
#ifdef CALIPER
  CALI_MARK_END("copy_particles_to_GPU");
#endif
#endif // CUDA

#ifdef COLLISION_DETECTION
  prepare_local_collision_queue();
#endif
  BondBreakage::clear_queue();
  auto particles = cell_structure.local_particles();
  auto ghost_particles = cell_structure.ghost_particles();
#ifdef ELECTROSTATICS
  if (system.coulomb.impl->extension) {
    if (auto icc = std::get_if<std::shared_ptr<ICCStar>>(
            get_ptr(system.coulomb.impl->extension))) {
      (**icc).iteration(cell_structure, particles, ghost_particles);
    }
  }
#endif
  init_forces(particles, ghost_particles, time_step, kT);

  calc_long_range_forces(particles);

  auto const elc_kernel = system.coulomb.pair_force_elc_kernel();
  auto const coulomb_kernel = system.coulomb.pair_force_kernel();
  auto const dipoles_kernel = system.dipoles.pair_force_kernel();

#ifdef ELECTROSTATICS
  auto const coulomb_cutoff = system.coulomb.cutoff();
#else
  auto const coulomb_cutoff = INACTIVE_CUTOFF;
#endif

#ifdef DIPOLES
  auto const dipole_cutoff = system.dipoles.cutoff();
#else
  auto const dipole_cutoff = INACTIVE_CUTOFF;
#endif

  short_range_loop(
      [coulomb_kernel_ptr = get_ptr(coulomb_kernel),
       &box_geo](Particle &p1, int bond_id, Utils::Span<Particle *> partners) {
        return add_bonded_force(p1, bond_id, partners, box_geo,
                                coulomb_kernel_ptr);
      },
      [coulomb_kernel_ptr = get_ptr(coulomb_kernel),
       dipoles_kernel_ptr = get_ptr(dipoles_kernel),
       elc_kernel_ptr = get_ptr(elc_kernel),
       &nonbonded_ias](Particle &p1, Particle &p2, Distance const &d) {
        auto const &ia_params =
            nonbonded_ias.get_ia_param(p1.type(), p2.type());
        add_non_bonded_pair_force(p1, p2, d.vec21, sqrt(d.dist2), d.dist2,
                                  ia_params, coulomb_kernel_ptr,
                                  dipoles_kernel_ptr, elc_kernel_ptr);
#ifdef COLLISION_DETECTION
        if (collision_params.mode != CollisionModeType::OFF)
          detect_collision(p1, p2, d.dist2);
#endif
      },
      cell_structure, maximal_cutoff(), maximal_cutoff_bonded(),
      VerletCriterion<>{system, skin, interaction_range(), coulomb_cutoff,
                        dipole_cutoff, collision_detection_cutoff()});

  Constraints::constraints.add_forces(box_geo, particles, get_sim_time());

  if (max_oif_objects) {
    // There are two global quantities that need to be evaluated:
    // object's surface and object's volume.
    for (int i = 0; i < max_oif_objects; i++) {
      auto const area_volume = boost::mpi::all_reduce(
          comm_cart, calc_oif_global(i, box_geo, cell_structure),
          std::plus<>());
      auto const oif_part_area = std::abs(area_volume[0]);
      auto const oif_part_vol = std::abs(area_volume[1]);
      if (oif_part_area < 1e-100 and oif_part_vol < 1e-100) {
        break;
      }
      add_oif_global_forces(area_volume, i, box_geo, cell_structure);
    }
  }

  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries.volume_conservation(cell_structure);

  if (system.lb.is_solver_set()) {
    LB::couple_particles(thermo_virtual, particles, ghost_particles, time_step);
  }

#ifdef CUDA
#ifdef CALIPER
  CALI_MARK_BEGIN("copy_forces_from_GPU");
#endif
  system.gpu.copy_forces_to_host(particles, this_node);
#ifdef CALIPER
  CALI_MARK_END("copy_forces_from_GPU");
#endif
#endif // CUDA

// VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
  virtual_sites()->back_transfer_forces_and_torques();
#endif

  // Communication Step: ghost forces
  cell_structure.ghosts_reduce_forces();

  // should be pretty late, since it needs to zero out the total force
  comfixed->apply(comm_cart, particles);

  // Needs to be the last one to be effective
  forcecap_cap(particles);

  // mark that forces are now up-to-date
  recalc_forces = false;
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
