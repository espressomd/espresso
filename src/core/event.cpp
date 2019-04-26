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
 *  Hook procedures.
 *
 *  Implemetation of event.hpp.
 */
#include "event.hpp"

#include "bonded_interactions/thermalized_bond.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "dpd.hpp"
#include "energy.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "ghosts.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "metadynamics.hpp"
#include "npt.hpp"
#include "nsquare.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "random.hpp"
#include "rattle.hpp"
#include "reaction_ensemble.hpp"
#include "rotation.hpp"
#include "statistics.hpp"
#include "swimmer_reaction.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include <utils/mpi/all_compare.hpp>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#ifdef SCAFACOS
#include "electrostatics_magnetostatics/scafacos.hpp"
#endif

/** whether the thermostat has to be reinitialized before integration */
static int reinit_thermo = 1;
static int reinit_electrostatics = 0;
static int reinit_magnetostatics = 0;

#ifdef CUDA
static int reinit_particle_comm_gpu = 1;
#endif

void on_program_start() {
  EVENT_TRACE(fprintf(stderr, "%d: on_program_start\n", this_node));

#ifdef CUDA
  cuda_init();
#endif

  /*
    call the initialization of the modules here
  */
  Random::init_random();

  init_node_grid();

  /* initially go for domain decomposition */
  topology_init(CELL_STRUCTURE_DOMDEC, &local_cells);

#ifdef LB_GPU
  if (this_node == 0) {
    //   lb_pre_init_gpu();
  }
#endif

#ifdef SWIMMER_REACTIONS
  reaction.eq_rate = 0.0;
  reaction.sing_mult = 0;
  reaction.swap = 0;
#endif

  /*
    call all initializations to do only on the master node here.
  */
  if (this_node == 0) {
    /* interaction_data.c: make sure 0<->0 ia always exists */
    make_particle_type_exist(0);
  }
}

void on_integration_start() {
  EVENT_TRACE(fprintf(stderr, "%d: on_integration_start\n", this_node));
  INTEG_TRACE(fprintf(
      stderr,
      "%d: on_integration_start: reinit_thermo = %d, resort_particles=%d\n",
      this_node, reinit_thermo, get_resort_particles()));

  /********************************************/
  /* sanity checks                            */
  /********************************************/

  integrator_sanity_checks();
#ifdef NPT
  integrator_npt_sanity_checks();
#endif
  interactions_sanity_checks();
#ifdef SWIMMER_REACTIONS
  reactions_sanity_checks();
#endif
#if defined(LB) || defined(LB_GPU)
  lb_lbfluid_on_integration_start();
#endif

  /********************************************/
  /* end sanity checks                        */
  /********************************************/
#ifdef CUDA
  if (reinit_particle_comm_gpu) {
    gpu_change_number_of_part_to_comm();
    reinit_particle_comm_gpu = 0;
  }
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
#endif

#ifdef METADYNAMICS
  meta_init();
#endif

  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermo_init();
    reinit_thermo = 0;
    recalc_forces = 1;
  }

  /* Ensemble preparation: NVT or NPT */
  integrate_ensemble_init();

  /* Update particle and observable information for routines in statistics.cpp
   */
  invalidate_obs();
  partCfg().invalidate();
  invalidate_fetch_cache();

#ifdef ADDITIONAL_CHECKS

  if (!Utils::Mpi::all_compare(comm_cart, cell_structure.type)) {
    runtimeErrorMsg() << "Nodes disagree about cell system type.";
  }

  if (!Utils::Mpi::all_compare(comm_cart, get_resort_particles())) {
    runtimeErrorMsg() << "Nodes disagree about resort type.";
  }

  if (!Utils::Mpi::all_compare(comm_cart, cell_structure.use_verlet_list)) {
    runtimeErrorMsg() << "Nodes disagree about use of verlet lists.";
  }

#ifdef ELECTROSTATICS
  if (!Utils::Mpi::all_compare(comm_cart, coulomb.method))
    runtimeErrorMsg() << "Nodes disagree about Coulomb long range method";
#endif
#ifdef DIPOLES
  if (!Utils::Mpi::all_compare(comm_cart, dipole.method))
    runtimeErrorMsg() << "Nodes disagree about dipolar long range method";
#endif
  check_global_consistency();
#endif /* ADDITIONAL_CHECKS */

  on_observable_calc();
}

void on_observable_calc() {
  EVENT_TRACE(fprintf(stderr, "%d: on_observable_calc\n", this_node));
  /* Prepare particle structure: Communication step: number of ghosts and ghost
   * information */

  cells_update_ghosts();

#ifdef ELECTROSTATICS
  if (reinit_electrostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_electrostatics\n", this_node));

    Coulomb::on_observable_calc();

    reinit_electrostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  if (reinit_magnetostatics) {
    EVENT_TRACE(fprintf(stderr, "%d: reinit_magnetostatics\n", this_node));

    Dipole::on_observable_calc();

    reinit_magnetostatics = 0;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef ELECTROKINETICS
  if (ek_initialized) {
    ek_integrate_electrostatics();
  }
#endif
}

void on_particle_charge_change() {
  reinit_electrostatics = 1;
  invalidate_obs();

  /* the particle information is no longer valid */
  partCfg().invalidate();
}

void on_particle_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_particle_change\n", this_node));

  set_resort_particles(Cells::RESORT_LOCAL);
  reinit_electrostatics = 1;
  reinit_magnetostatics = 1;

#ifdef LB_GPU
  lb_lbfluid_invalidate_particle_allocation();
#endif
#ifdef CUDA
  reinit_particle_comm_gpu = 1;
#endif
  invalidate_obs();

  /* the particle information is no longer valid */
  partCfg().invalidate();

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
}

void on_coulomb_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_coulomb_change\n", this_node));
  invalidate_obs();

#ifdef ELECTROSTATICS
  Coulomb::on_coulomb_change();
#endif /* ELECTROSTATICS */

#ifdef DIPOLES
  Dipole::on_coulomb_change();
#endif /* ifdef DIPOLES */

  /* all Coulomb methods have a short range part, aka near field
     correction. Even in case of switching off, we should call this,
     since the required cutoff might have reduced. */
  on_short_range_ia_change();

#ifdef CUDA
  reinit_particle_comm_gpu = 1;
#endif
  recalc_forces = 1;
}

void on_short_range_ia_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_short_range_ia_changes\n", this_node));
  invalidate_obs();

  recalc_maximal_cutoff();
  cells_on_geometry_change(0);

  recalc_forces = 1;
}

void on_constraint_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_constraint_change\n", this_node));
  invalidate_obs();
  recalc_forces = 1;
}

void on_lbboundary_change() {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  EVENT_TRACE(fprintf(stderr, "%d: on_lbboundary_change\n", this_node));
  invalidate_obs();

  LBBoundaries::lb_init_boundaries();

  recalc_forces = 1;
#endif
}

void on_resort_particles() {
  EVENT_TRACE(fprintf(stderr, "%d: on_resort_particles\n", this_node));
#ifdef ELECTROSTATICS
  Coulomb::on_resort_particles();
#endif /* ifdef ELECTROSTATICS */

  /* DIPOLAR interactions so far don't need this */

  recalc_forces = 1;
}

void on_boxl_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_boxl_change\n", this_node));

  grid_changed_box_l();
  /* Electrostatics cutoffs mostly depend on the system size,
     therefore recalculate them. */
  recalc_maximal_cutoff();
  cells_on_geometry_change(0);

/* Now give methods a chance to react to the change in box length */
#ifdef ELECTROSTATICS
  Coulomb::on_boxl_change();
#endif

#ifdef DIPOLES
  Dipole::on_boxl_change();
#endif

#ifdef LB
  lb_lbfluid_init();
#ifdef LB_BOUNDARIES
  LBBoundaries::lb_init_boundaries();
#endif
#endif
}

void on_cell_structure_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_cell_structure_change\n", this_node));

/* Now give methods a chance to react to the change in cell
   structure.  Most ES methods need to reinitialize, as they depend
   on skin, node grid and so on. Only for a change in box length we
   have separate, faster methods, as this might happen frequently
   in a NpT simulation. */
#ifdef ELECTROSTATICS
  Coulomb::init();
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  Dipole::init();
#endif /* ifdef DIPOLES */

#ifdef LB
  if (lattice_switch == ActiveLB::CPU) {
    runtimeErrorMsg()
        << "The CPU LB does not currently support handling changes of the MD "
           "cell geometry. Setup the cell system, skin and interactions before "
           "activating the CPU LB.";
  }
#endif
}

void on_temperature_change() {
#if defined(LB) || defined(LB_GPU)
  lb_lbfluid_reinit_parameters();
#endif
}

void on_parameter_change(int field) {
  EVENT_TRACE(
      fprintf(stderr, "%d: shon_parameter_change %d\n", this_node, field));

  switch (field) {
  case FIELD_BOXL:
    on_boxl_change();
    break;
  case FIELD_MIN_GLOBAL_CUT:
    recalc_maximal_cutoff();
    cells_on_geometry_change(0);
    break;
  case FIELD_SKIN:
    cells_on_geometry_change(0);
    break;
  case FIELD_PERIODIC:
#ifdef SCAFACOS
#ifdef ELECTROSTATICS
    if (coulomb.method == COULOMB_SCAFACOS) {
      Scafacos::update_system_params();
    }
#endif
#ifdef DIPOLES
    if (dipole.method == DIPOLAR_SCAFACOS) {
      Scafacos::update_system_params();
    }
#endif

#endif
    cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
    break;
  case FIELD_NODEGRID:
    grid_changed_n_nodes();
    cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
    break;
  case FIELD_MINNUMCELLS:
  case FIELD_MAXNUMCELLS:
    cells_re_init(CELL_STRUCTURE_CURRENT);
    break;
  case FIELD_TEMPERATURE:
    on_temperature_change();
    reinit_thermo = 1;
    break;
  case FIELD_TIMESTEP:
#if defined(LB) || defined(LB_GPU)
    lb_lbfluid_reinit_parameters();
#endif
  case FIELD_LANGEVIN_GAMMA:
  case FIELD_LANGEVIN_GAMMA_ROTATION:
  case FIELD_NPTISO_G0:
  case FIELD_NPTISO_GV:
  case FIELD_NPTISO_PISTON:
    reinit_thermo = 1;
    break;
#ifdef NPT
  case FIELD_INTEG_SWITCH:
    if (integ_switch != INTEG_METHOD_NPT_ISO)
      nptiso.invalidate_p_vel = 1;
    break;
#endif
  case FIELD_THERMO_SWITCH:
    /* DPD needs ghost velocities, other thermostats not */
    on_ghost_flags_change();
#ifdef DPD
    if (not(thermo_switch & THERMO_DPD)) {
      dpd_switch_off();
    }
#endif
    break;
#ifdef LB
  case FIELD_LATTICE_SWITCH:
    /* LB needs ghost velocities */
    on_ghost_flags_change();
    break;
#endif
  case FIELD_FORCE_CAP:
    /* If the force cap changed, forces are invalid */
    invalidate_obs();
    recalc_forces = 1;
    break;
  case FIELD_RIGIDBONDS:
    /* Rattle bonds needs ghost velocities */
    on_ghost_flags_change();
    break;
  case FIELD_THERMALIZEDBONDS:
    /* Thermalized distance bonds needs ghost velocities */
    on_ghost_flags_change();
    break;
  case FIELD_SIMTIME:
    recalc_forces = 1;
    break;
  }
}

void on_ghost_flags_change() {
  EVENT_TRACE(fprintf(stderr, "%d: on_ghost_flags_change\n", this_node));
  /* that's all we change here */
  extern int ghosts_have_v;
  extern int ghosts_have_bonds;

  ghosts_have_v = 0;
  ghosts_have_bonds = 0;

/* DPD and LB need also ghost velocities */
#ifdef LB
  if (lattice_switch == ActiveLB::CPU)
    ghosts_have_v = 1;
#endif
#ifdef BOND_CONSTRAINT
  if (n_rigidbonds)
    ghosts_have_v = 1;
#endif
#ifdef DPD
  // maybe we have to add a new global to differ between compile in and actual
  // use.
  if (thermo_switch & THERMO_DPD)
    ghosts_have_v = 1;
#endif
#ifdef VIRTUAL_SITES
  // If they have velocities, VIRUTAL_SITES need v to update v of virtual sites
  if (virtual_sites()->get_have_velocity()) {
    ghosts_have_v = 1;
  };
#endif
  // THERMALIZED_DIST_BOND needs v to calculate v_com and v_dist for thermostats
  if (n_thermalized_bonds) {
    ghosts_have_v = 1;
    ghosts_have_bonds = 1;
  }
#ifdef COLLISION_DETECTION
  if (collision_params.mode) {
    ghosts_have_bonds = 1;
  }
#endif
}
