
/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *  Force calculation.
 *
 *  The corresponding header file is forces.hpp.
 */

#include "EspressoSystemInterface.hpp"

#include "comfixed_global.hpp"
#include "constraints.hpp"
#include "electrostatics_magnetostatics/icc.hpp"
#include "electrostatics_magnetostatics/p3m_gpu.hpp"
#include "forcecap.hpp"
#include "forces_inline.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "immersed_boundaries.hpp"
#include "short_range_loop.hpp"

#include <profiler/profiler.hpp>

#include <cassert>

ActorList forceActors;

void init_forces() {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  /* The force initialization depends on the used thermostat and the
     thermodynamic ensemble */

#ifdef NPT
  /* reset virial part of instantaneous pressure */
  if (integ_switch == INTEG_METHOD_NPT_ISO)
    nptiso.p_vir[0] = nptiso.p_vir[1] = nptiso.p_vir[2] = 0.0;
#endif

  /* initialize forces with Langevin thermostat forces
     or zero depending on the thermostat
     set torque to zero for all and rescale quaternions
  */
  for (auto &p : local_cells.particles()) {
    init_local_particle_force(&p);
  }

  /* initialize ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  for (auto &p : ghost_cells.particles()) {
    init_ghost_force(&p);
  }
}

void init_forces_ghosts() {
  for (auto &p : ghost_cells.particles()) {
    init_ghost_force(&p);
  }
}

// This function is no longer called from force_calc().
// The check was moved to rescale_fores() to avoid an additional iteration over
// all particles
void check_forces() {
  for (auto &p : local_cells.particles()) {
    check_particle_force(&p);
  }

  for (auto &p : ghost_cells.particles()) {
    check_particle_force(&p);
  }
}

void force_calc() {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  espressoSystemInterface.update();

#ifdef COLLISION_DETECTION
  prepare_local_collision_queue();
#endif

#ifdef ELECTROSTATICS
  iccp3m_iteration();
#endif
  init_forces();

  for (auto actor = forceActors.begin(); actor != forceActors.end(); ++actor) {
    (*actor)->computeForces(espressoSystemInterface);
#ifdef ROTATION
    (*actor)->computeTorques(espressoSystemInterface);
#endif
  }

  calc_long_range_forces();

  // Only calculate pair forces if the maximum cutoff is >0
  if (max_cut > 0) {
    short_range_loop([](Particle &p) { add_single_particle_force(&p); },
                     [](Particle &p1, Particle &p2, Distance &d) {
                       add_non_bonded_pair_force(&(p1), &(p2), d.vec21.data(),
                                                 sqrt(d.dist2), d.dist2);
                     });
  } else {
    // Otherwise only do single-particle contributions
    for (auto &p : local_cells.particles()) {
      add_single_particle_force(&p);
    }
  }
  auto local_parts = local_cells.particles();
  Constraints::constraints.add_forces(local_parts, sim_time);

#ifdef OIF_GLOBAL_FORCES
  if (max_oif_objects) {
    double area_volume[2]; // There are two global quantities that need to be
    // evaluated: object's surface and object's volume. One
    // can add another quantity.
    area_volume[0] = 0.0;
    area_volume[1] = 0.0;
    for (int i = 0; i < max_oif_objects; i++) {
      calc_oif_global(area_volume, i);
      if (fabs(area_volume[0]) < 1e-100 && fabs(area_volume[1]) < 1e-100)
        break;
      add_oif_global_forces(area_volume, i);
    }
  }
#endif

#ifdef IMMERSED_BOUNDARY
  // Must be done here. Forces need to be ghost-communicated
  immersed_boundaries.volume_conservation();
#endif

#if defined(LB_GPU) || defined(LB)
  lb_lbcoupling_calc_particle_lattice_ia(thermo_virtual);
#endif

#ifdef METADYNAMICS
  /* Metadynamics main function */
  meta_perform();
#endif

#ifdef CUDA
  copy_forces_from_GPU(local_cells.particles());
#endif

// VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
  if (virtual_sites()->need_ghost_comm_before_back_transfer()) {
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    init_forces_ghosts();
  }
  virtual_sites()->back_transfer_forces_and_torques();
#endif

  // Communication Step: ghost forces
  ghost_communicator(&cell_structure.collect_ghost_force_comm);

  auto local_particles = local_cells.particles();
  // should be pretty late, since it needs to zero out the total force
  comfixed.apply(comm_cart, local_particles);

  // Needs to be the last one to be effective
  forcecap_cap(local_cells.particles());

  // mark that forces are now up-to-date
  recalc_forces = 0;
}

void calc_long_range_forces() {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (elc_params.dielectric_contrast_on) {
      ELC_P3M_modify_p3m_sums_both();
      ELC_p3m_charge_assign_both();
      ELC_P3M_self_forces();
    } else
      p3m_charge_assign();

    p3m_calc_kspace_forces(1, 0);

    if (elc_params.dielectric_contrast_on)
      ELC_P3M_restore_p3m_sums();

    ELC_add_force();

    break;
#endif
#ifdef CUDA
  case COULOMB_P3M_GPU:
    if (this_node == 0) {
      FORCE_TRACE(printf("Computing GPU P3M forces.\n"));
      p3m_gpu_add_farfield_force();
    }
    /* there is no NPT handling here as long as we cannot compute energies.
       This is checked in integrator_npt_sanity_checks() when integration
       starts. */
    break;
#endif
#ifdef P3M
  case COULOMB_P3M:
    FORCE_TRACE(printf("%d: Computing P3M forces.\n", this_node));
    p3m_charge_assign();
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += p3m_calc_kspace_forces(1, 1);
    else
#endif
      p3m_calc_kspace_forces(1, 0);
    break;
#endif
  case COULOMB_MMM2D:
    MMM2D_add_far_force();
    MMM2D_dielectric_layers_force_contribution();
    break;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    assert(!Scafacos::dipolar());
    Scafacos::add_long_range_force();
    break;
#endif
  default:
    break;
  }

#ifdef ELECTROKINETICS
  /* If enabled, calculate electrostatics contribution from electrokinetics
   * species. */
  ek_calculate_electrostatic_coupling();
#endif

#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* calculate k-space part of the magnetostatic interaction. */
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    add_mdlc_force_corrections();
  // fall through
  case DIPOLAR_P3M:
    dp3m_dipole_assign();
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      nptiso.p_vir[0] += dp3m_calc_kspace_forces(1, 1);
      fprintf(stderr, "dipolar_P3M at this moment is added to p_vir[0]\n");
    } else
#endif
      dp3m_calc_kspace_forces(1, 0);

    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    dawaanr_calculations(1, 0);
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_DS:
    add_mdlc_force_corrections();
    // fall through
#endif
  case DIPOLAR_DS:
    magnetic_dipolar_direct_sum_calculations(1, 0);
    break;
  case DIPOLAR_DS_GPU:
    // Do nothing. It's an actor
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    // Do nothing, it's an actor.
    break;
#endif // BARNES_HUT
#ifdef SCAFACOS_DIPOLES
  case DIPOLAR_SCAFACOS:
    assert(Scafacos::dipolar());
    Scafacos::add_long_range_force();
#endif
  case DIPOLAR_NONE:
    break;
  default:
    runtimeErrorMsg() << "unknown dipolar method";
    break;
  }
#endif /*ifdef DIPOLES */
}
