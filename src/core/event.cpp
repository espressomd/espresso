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
 *  Hook procedures.
 *
 *  Implementation of event.hpp.
 */
#include "event.hpp"

#include "Particle.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "immersed_boundaries.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "statistics.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include <utils/mpi/all_compare.hpp>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#ifdef SCAFACOS
#include "electrostatics_magnetostatics/scafacos.hpp"
#endif

/** whether the thermostat has to be reinitialized before integration */
static bool reinit_thermo = true;
static int reinit_electrostatics = false;
static int reinit_magnetostatics = false;

#if defined(OPEN_MPI) &&                                                       \
    (OMPI_MAJOR_VERSION == 2 && OMPI_MINOR_VERSION <= 1 ||                     \
     OMPI_MAJOR_VERSION == 3 &&                                                \
         (OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION <= 2 ||              \
          OMPI_MINOR_VERSION == 1 && OMPI_RELEASE_VERSION <= 2))
/** Workaround for segmentation fault "Signal code: Address not mapped (1)"
 *  that happens when the visualizer is used. This is a bug in OpenMPI 2.0-2.1,
 *  3.0.0-3.0.2 and 3.1.0-3.1.2
 *  https://github.com/espressomd/espresso/issues/3056
 */
#define OPENMPI_BUG_MPI_ALLOC_MEM
#endif

void on_program_start() {
#ifdef CUDA
  cuda_init();
#endif

  init_node_grid();

  /* initially go for domain decomposition */
  cells_re_init(CELL_STRUCTURE_DOMDEC, INACTIVE_CUTOFF);

  /*
    call all initializations to do only on the master node here.
  */
  if (this_node == 0) {
    /* interaction_data.c: make sure 0<->0 ia always exists */
    make_particle_type_exist(0);
  }
}

void on_integration_start() {
  /********************************************/
  /* sanity checks                            */
  /********************************************/

  integrator_sanity_checks();
#ifdef NPT
  integrator_npt_sanity_checks();
#endif
  interactions_sanity_checks();
  lb_lbfluid_on_integration_start();

  /********************************************/
  /* end sanity checks                        */
  /********************************************/

#ifdef CUDA
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
#endif

  // Here we initialize volume conservation
  // This function checks if the reference volumes have been set and if
  // necessary calculates them
  immersed_boundaries.init_volume_conservation(cell_structure);

  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermo_init();
    reinit_thermo = false;
    recalc_forces = true;
  }

#ifdef NPT
  npt_ensemble_init(box_geo);
#endif

  invalidate_obs();
  partCfg().invalidate();
  invalidate_fetch_cache();

#ifdef ADDITIONAL_CHECKS

  if (!Utils::Mpi::all_compare(comm_cart, cell_structure.type)) {
    runtimeErrorMsg() << "Nodes disagree about cell system type.";
  }

  if (!Utils::Mpi::all_compare(comm_cart, cell_structure.use_verlet_list)) {
    runtimeErrorMsg() << "Nodes disagree about use of verlet lists.";
  }
#ifndef OPENMPI_BUG_MPI_ALLOC_MEM
#ifdef ELECTROSTATICS
  if (!Utils::Mpi::all_compare(comm_cart, coulomb.method))
    runtimeErrorMsg() << "Nodes disagree about Coulomb long range method";
#endif
#ifdef DIPOLES
  if (!Utils::Mpi::all_compare(comm_cart, dipole.method))
    runtimeErrorMsg() << "Nodes disagree about dipolar long range method";
#endif
#endif
  check_global_consistency();
#endif /* ADDITIONAL_CHECKS */

  on_observable_calc();
}

void on_observable_calc() {
  /* Prepare particle structure: Communication step: number of ghosts and ghost
   * information */
  cells_update_ghosts(global_ghost_flags());
  update_dependent_particles();
#ifdef ELECTROSTATICS
  if (reinit_electrostatics) {
    Coulomb::on_observable_calc();
    reinit_electrostatics = false;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  if (reinit_magnetostatics) {
    Dipole::on_observable_calc();
    reinit_magnetostatics = false;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef ELECTROKINETICS
  if (ek_initialized) {
    ek_integrate_electrostatics();
  }
#endif

  clear_particle_node();
}

void on_particle_charge_change() {
  reinit_electrostatics = true;
  invalidate_obs();

  /* the particle information is no longer valid */
  partCfg().invalidate();
}

void on_particle_change() {
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  reinit_electrostatics = true;
  reinit_magnetostatics = true;

  invalidate_obs();

  /* the particle information is no longer valid */
  partCfg().invalidate();

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
}

void on_coulomb_change() {
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

  recalc_forces = true;
}

void on_short_range_ia_change() {
  invalidate_obs();

  cells_on_geometry_change(false);

  recalc_forces = true;
}

void on_constraint_change() {
  invalidate_obs();
  recalc_forces = true;
}

void on_lbboundary_change() {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  invalidate_obs();

  LBBoundaries::lb_init_boundaries();

  recalc_forces = true;
#endif
}

void on_resort_particles() { recalc_forces = true; }

void on_boxl_change() {
  grid_changed_box_l(box_geo);
  /* Electrostatics cutoffs mostly depend on the system size,
     therefore recalculate them. */
  cells_on_geometry_change(false);

/* Now give methods a chance to react to the change in box length */
#ifdef ELECTROSTATICS
  Coulomb::on_boxl_change();
#endif

#ifdef DIPOLES
  Dipole::on_boxl_change();
#endif

  lb_lbfluid_init();
#ifdef LB_BOUNDARIES
  LBBoundaries::lb_init_boundaries();
#endif
}

void on_cell_structure_change() {
/* Now give methods a chance to react to the change in cell
   structure. Most ES methods need to reinitialize, as they depend
   on skin, node grid and so on. Only for a change in box length we
   have separate, faster methods, as this might happen frequently
   in a NpT simulation. */
#ifdef ELECTROSTATICS
  Coulomb::init();
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  Dipole::init();
#endif /* ifdef DIPOLES */

  if (lattice_switch == ActiveLB::CPU) {
    runtimeErrorMsg()
        << "The CPU LB does not currently support handling changes of the MD "
           "cell geometry. Setup the cell system, skin and interactions before "
           "activating the CPU LB.";
  }
}

void on_temperature_change() { lb_lbfluid_reinit_parameters(); }

void on_parameter_change(int field) {

  switch (field) {
  case FIELD_BOXL:
    on_boxl_change();
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
  case FIELD_MIN_GLOBAL_CUT:
  case FIELD_SKIN:
    cells_on_geometry_change(false);
    break;
  case FIELD_NODEGRID:
    grid_changed_n_nodes();
    cells_re_init(CELL_STRUCTURE_CURRENT, cell_structure.min_range);
    break;
  case FIELD_TEMPERATURE:
    on_temperature_change();
    reinit_thermo = true;
    break;
  case FIELD_TIMESTEP:
    lb_lbfluid_reinit_parameters();
  case FIELD_LANGEVIN_GAMMA:
  case FIELD_LANGEVIN_GAMMA_ROTATION:
  case FIELD_NPTISO_G0:
  case FIELD_NPTISO_GV:
  case FIELD_NPTISO_PISTON:
    reinit_thermo = true;
    break;
#ifdef NPT
  case FIELD_INTEG_SWITCH:
    if (integ_switch != INTEG_METHOD_NPT_ISO)
      nptiso.invalidate_p_vel = true;
    break;
#endif
    break;
  case FIELD_FORCE_CAP:
    /* If the force cap changed, forces are invalid */
    invalidate_obs();
    recalc_forces = true;
    break;
  case FIELD_THERMO_SWITCH:
  case FIELD_LATTICE_SWITCH:
  case FIELD_RIGIDBONDS:
  case FIELD_THERMALIZEDBONDS:
    break;
  case FIELD_SIMTIME:
    recalc_forces = true;
    break;
  }
}

/**
 * @brief Returns the ghost flags required for running pair
 *        kernels for the global state, e.g. the force calculation.
 * @return Required data parts;
 */
unsigned global_ghost_flags() {
  /* Position and Properties are always requested. */
  unsigned data_parts = Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES;

  if (lattice_switch == ActiveLB::CPU)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (thermo_switch & THERMO_DPD)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (n_thermalized_bonds) {
    data_parts |= Cells::DATA_PART_MOMENTUM;
    data_parts |= Cells::DATA_PART_BONDS;
  }

#ifdef COLLISION_DETECTION
  if (collision_params.mode) {
    data_parts |= Cells::DATA_PART_BONDS;
  }
#endif

  return data_parts;
}

void update_dependent_particles() {
#ifdef VIRTUAL_SITES
  virtual_sites()->update();
  cells_update_ghosts(global_ghost_flags());
#endif

#ifdef ELECTROSTATICS
  Coulomb::update_dependent_particles();
#endif
}
