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

#include "bonded_interactions/thermalized_bond.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "immersed_boundaries.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#ifdef SCAFACOS
#include "electrostatics_magnetostatics/scafacos.hpp"
#endif

#include <utils/mpi/all_compare.hpp>

#include <cstdio>

#include <mpi.h>

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
  if (this_node == 0) {
    try {
      cuda_init();
    } catch (cuda_runtime_error const &err) {
      // pass
    }
  }
#endif

  init_node_grid();

  /* initially go for domain decomposition */
  cells_re_init(CELL_STRUCTURE_DOMDEC);

  if (this_node == 0) {
    /* make sure interaction 0<->0 always exists */
    make_particle_type_exist(0);
  }
}

void on_integration_start(double time_step) {
  /********************************************/
  /* sanity checks                            */
  /********************************************/

  integrator_sanity_checks();
#ifdef NPT
  integrator_npt_sanity_checks();
#endif
  if (long_range_interactions_sanity_checks()) {
    runtimeErrorMsg() << "Long-range interactions returned an error.";
  }
  lb_lbfluid_sanity_checks(time_step);

  /********************************************/
  /* end sanity checks                        */
  /********************************************/

  lb_lbfluid_on_integration_start();

#ifdef CUDA
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
#endif

  /* Prepare the thermostat */
  if (reinit_thermo) {
    thermo_init(time_step);
    reinit_thermo = false;
    recalc_forces = true;
  }

#ifdef NPT
  npt_ensemble_init(box_geo);
#endif

  partCfg().invalidate();
  invalidate_fetch_cache();

#ifdef ADDITIONAL_CHECKS
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
#endif /* ELECTROSTATICS */

#ifdef DIPOLES
  if (reinit_magnetostatics) {
    Dipole::on_observable_calc();
    reinit_magnetostatics = false;
  }
#endif /* DIPOLES */

#ifdef ELECTROKINETICS
  if (ek_initialized) {
    ek_integrate_electrostatics();
  }
#endif /* ELECTROKINETICS */

  clear_particle_node();
}

void on_particle_charge_change() {
  reinit_electrostatics = true;

  /* the particle information is no longer valid */
  partCfg().invalidate();
}

void on_particle_change() {
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  reinit_electrostatics = true;
  reinit_magnetostatics = true;
  recalc_forces = true;

  /* the particle information is no longer valid */
  partCfg().invalidate();

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
}

void on_coulomb_change() {

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
  cells_re_init(cell_structure.decomposition_type());

  recalc_forces = true;
}

void on_constraint_change() { recalc_forces = true; }

void on_lbboundary_change() {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  LBBoundaries::lb_init_boundaries();

  recalc_forces = true;
#endif
}

void on_boxl_change(bool skip_method_adaption) {
  grid_changed_box_l(box_geo);
  /* Electrostatics cutoffs mostly depend on the system size,
   * therefore recalculate them. */
  cells_re_init(cell_structure.decomposition_type());

  if (not skip_method_adaption) {
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
}

void on_cell_structure_change() {
  clear_particle_node();

  /* Now give methods a chance to react to the change in cell structure.
   * Most ES methods need to reinitialize, as they depend on skin,
   * node grid and so on. */
#ifdef ELECTROSTATICS
  Coulomb::init();
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  Dipole::init();
#endif /* ifdef DIPOLES */
}

void on_temperature_change() { lb_lbfluid_reinit_parameters(); }

void on_periodicity_change() {
#ifdef SCAFACOS
#ifdef ELECTROSTATICS
  if (coulomb.method == COULOMB_SCAFACOS) {
    Scafacos::fcs_coulomb()->update_system_params();
  }
#endif
#ifdef SCAFACOS_DIPOLES
  if (dipole.method == DIPOLAR_SCAFACOS) {
    Scafacos::fcs_dipoles()->update_system_params();
  }
#endif
#endif
#ifdef STOKESIAN_DYNAMICS
  if (integ_switch == INTEG_METHOD_SD) {
    if (box_geo.periodic(0) || box_geo.periodic(1) || box_geo.periodic(2))
      runtimeErrorMsg() << "Illegal box periodicity for Stokesian Dynamics: "
                        << box_geo.periodic(0) << " " << box_geo.periodic(1)
                        << " " << box_geo.periodic(2) << "\n"
                        << "  Required: 0 0 0\n";
  }
#endif
  on_skin_change();
}

void on_skin_change() {
  cells_re_init(cell_structure.decomposition_type());
  on_coulomb_change();
}

void on_thermostat_param_change() { reinit_thermo = true; }

void on_timestep_change() {
  lb_lbfluid_reinit_parameters();
  on_thermostat_param_change();
}

void on_simtime_change() { recalc_forces = true; }

void on_forcecap_change() { recalc_forces = true; }

void on_nodegrid_change() {
  grid_changed_n_nodes();
  cells_re_init(cell_structure.decomposition_type());
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

  // Here we initialize volume conservation
  // This function checks if the reference volumes have been set and if
  // necessary calculates them
  immersed_boundaries.init_volume_conservation(cell_structure);
}
