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
 *  Hook procedures.
 *
 *  Implementation of event.hpp.
 */
#include "event.hpp"

#include "bonded_interactions/thermalized_bond.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "electrostatics/coulomb.hpp"
#include "electrostatics/icc.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "immersed_boundaries.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_node.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include <utils/mpi/all_compare.hpp>

#include <mpi.h>

/** whether the thermostat has to be reinitialized before integration */
static bool reinit_thermo = true;
#ifdef ELECTROSTATICS
/** whether electrostatics actor has to be reinitialized on observable calc */
static bool reinit_electrostatics = false;
#endif
#ifdef DIPOLES
/** whether magnetostatics actor has to be reinitialized on observable calc */
static bool reinit_magnetostatics = false;
#endif

void on_program_start() {
#ifdef CUDA
  if (this_node == 0) {
    try {
      cuda_init();
    } catch (cuda_runtime_error const &) {
      // pass
    }
  }
#endif

  init_node_grid();

  /* initially go for regular decomposition */
  cells_re_init(CellStructureType::CELL_STRUCTURE_REGULAR);

  /* make sure interaction 0<->0 always exists */
  make_particle_type_exist_local(0);
}

void on_integration_start(double time_step) {
  /********************************************/
  /* sanity checks                            */
  /********************************************/

  integrator_sanity_checks();
#ifdef NPT
  integrator_npt_sanity_checks();
#endif
  long_range_interactions_sanity_checks();
  LB::sanity_checks(time_step);

  /********************************************/
  /* end sanity checks                        */
  /********************************************/

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
#ifdef ELECTROSTATICS
  {
    auto const &actor = electrostatics_actor;
    if (!Utils::Mpi::all_compare(comm_cart, static_cast<bool>(actor)) or
        (actor and !Utils::Mpi::all_compare(comm_cart, (*actor).which())))
      runtimeErrorMsg() << "Nodes disagree about Coulomb long-range method";
  }
#endif
#ifdef DIPOLES
  {
    auto const &actor = magnetostatics_actor;
    if (!Utils::Mpi::all_compare(comm_cart, static_cast<bool>(actor)) or
        (actor and !Utils::Mpi::all_compare(comm_cart, (*actor).which())))
      runtimeErrorMsg() << "Nodes disagree about dipolar long-range method";
  }
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
    Dipoles::on_observable_calc();
    reinit_magnetostatics = false;
  }
#endif /* DIPOLES */

  clear_particle_node();
}

void on_particle_charge_change() {
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
#endif

  /* the particle information is no longer valid */
  partCfg().invalidate();
}

void on_particle_change() {
  if (cell_structure.decomposition_type() ==
      CellStructureType::CELL_STRUCTURE_HYBRID) {
    cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  } else {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
#endif
#ifdef DIPOLES
  reinit_magnetostatics = true;
#endif
  recalc_forces = true;

  /* the particle information is no longer valid */
  partCfg().invalidate();

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
}

void on_coulomb_and_dipoles_change() {
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
  Coulomb::on_coulomb_change();
#endif
#ifdef DIPOLES
  reinit_magnetostatics = true;
  Dipoles::on_dipoles_change();
#endif
  on_short_range_ia_change();
}

void on_coulomb_change() {
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
  Coulomb::on_coulomb_change();
#endif
  on_short_range_ia_change();
}

void on_dipoles_change() {
#ifdef DIPOLES
  reinit_magnetostatics = true;
  Dipoles::on_dipoles_change();
#endif
  on_short_range_ia_change();
}

void on_non_bonded_ia_change() {
  maximal_cutoff_nonbonded();
  on_short_range_ia_change();
}

void on_short_range_ia_change() {
  cells_re_init(cell_structure.decomposition_type());
  recalc_forces = true;
}

void on_constraint_change() { recalc_forces = true; }

void on_lb_boundary_conditions_change() { recalc_forces = true; }

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
    Dipoles::on_boxl_change();
#endif

    LB::init();
  }
}

void on_cell_structure_change() {
  clear_particle_node();

  if (lattice_switch == ActiveLB::WALBERLA_LB) {
    throw std::runtime_error(
        "LB does not currently support handling changes of the MD cell "
        "geometry. Setup the cell system, skin and interactions before "
        "activating the CPU LB.");
  }

  /* Now give methods a chance to react to the change in cell structure.
   * Most ES methods need to reinitialize, as they depend on skin,
   * node grid and so on. */
#ifdef ELECTROSTATICS
  Coulomb::on_cell_structure_change();
#endif

#ifdef DIPOLES
  Dipoles::on_cell_structure_change();
#endif
}

void on_temperature_change() {
  if (lattice_switch != ActiveLB::NONE) {
    throw std::runtime_error("Temperature change not supported by LB");
  }
}

void on_periodicity_change() {
#ifdef ELECTROSTATICS
  Coulomb::on_periodicity_change();
#endif

#ifdef DIPOLES
  Dipoles::on_periodicity_change();
#endif

#ifdef STOKESIAN_DYNAMICS
  if (integ_switch == INTEG_METHOD_SD) {
    if (box_geo.periodic(0) || box_geo.periodic(1) || box_geo.periodic(2))
      runtimeErrorMsg() << "Stokesian Dynamics requires periodicity "
                        << "(False, False, False)\n";
  }
#endif
  on_skin_change();
}

void on_skin_change() {
  cells_re_init(cell_structure.decomposition_type());
  on_coulomb_and_dipoles_change();
}

void on_thermostat_param_change() { reinit_thermo = true; }

void on_timestep_change() {
  if (lattice_switch != ActiveLB::NONE) {
    throw std::runtime_error("Time step change not supported by LB");
  }
  on_thermostat_param_change();
}

void on_forcecap_change() { recalc_forces = true; }

void on_node_grid_change() {
  grid_changed_n_nodes();
#ifdef ELECTROSTATICS
  Coulomb::on_node_grid_change();
#endif
#ifdef DIPOLES
  Dipoles::on_node_grid_change();
#endif
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

  if (lattice_switch == ActiveLB::WALBERLA_LB)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (thermo_switch & THERMO_DPD)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (n_thermalized_bonds) {
    data_parts |= Cells::DATA_PART_MOMENTUM;
    data_parts |= Cells::DATA_PART_BONDS;
  }

#ifdef COLLISION_DETECTION
  if (collision_params.mode != CollisionModeType::OFF) {
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
  update_icc_particles();
#endif

  // Here we initialize volume conservation
  // This function checks if the reference volumes have been set and if
  // necessary calculates them
  immersed_boundaries.init_volume_conservation(cell_structure);
}
