/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _GLOBAL_HPP
#define _GLOBAL_HPP

/** \file global.hpp This file contains the code for access to globally
    defined variables using the script command setmd. Please refer to
    the Developer's guide, section "Adding global variables", for
    details on how to add new variables in the interpreter's
    space.  */

/**********************************************
 * description of global variables
 * add any variable that should be handled
 * automatically in global.cpp. This includes
 * distribution to other nodes and
 * read/user-defined access from Tcl.
 **********************************************/

/** Issue REQ_BCAST_PAR: broadcast a parameter from datafield.
         @param i the number from \ref global.hpp "global.hpp" referencing the
         datafield.
         @return nonzero on error
     */
int mpi_bcast_parameter(int i);

/*
 * @brief Check if all the global fields a synchronized between
 *        the nodes.
 */
void check_global_consistency();

/** \name Field Enumeration
    These numbers identify the variables given in
    \ref #fields for use with \ref mpi_bcast_parameter.
*/
/*@{*/
enum Fields {
  /** index of \ref box_l in \ref #fields */
  FIELD_BOXL = 0,
  /** index of \ref DomainDecomposition::cell_grid in  \ref #fields */
  FIELD_CELLGRID,
  /** index of \ref DomainDecomposition::cell_size in  \ref #fields */
  FIELD_CELLSIZE,
  /** index of \ref langevin_gamma in  \ref #fields */
  FIELD_LANGEVIN_GAMMA,
  /** index of \ref lees_edwards_offset in \ref #fields */
  FIELD_LEES_EDWARDS_OFFSET,
  /** index of \ref integ_switch in \ref #fields */
  FIELD_INTEG_SWITCH,
  /** index of \ref local_box_l in \ref #fields */
  FIELD_LBOXL,
  /** index of \ref max_cut in \ref #fields */
  FIELD_MCUT,
  /** index of \ref max_num_cells  in \ref #fields */
  FIELD_MAXNUMCELLS,
  /** index of \ref max_seen_particle in \ref #fields */
  FIELD_MAXPART,
  /** index of \ref max_range in \ref #fields */
  FIELD_MAXRANGE,
  /** index of \ref max_skin in  \ref #fields */
  FIELD_MAXSKIN,
  /** index of \ref min_num_cells  in \ref #fields */
  FIELD_MINNUMCELLS,
  /** index of \ref n_layers in  \ref #fields */
  FIELD_NLAYERS,
  /** index of \ref n_nodes in \ref #fields */
  FIELD_NNODES,
  /** index of \ref n_part in  \ref #fields */
  FIELD_NPART,
  /** index of \ref n_particle_types in \ref #fields */
  FIELD_NPARTTYPE,
  /** index of \ref n_rigidbonds in \ref #fields */
  FIELD_RIGIDBONDS,
  /** index of \ref node_grid in \ref #fields */
  FIELD_NODEGRID,
  /** index of \ref nptiso_gamma in \ref #fields */
  FIELD_NPTISO_G0,
  /** index of \ref nptiso_gammav in \ref #fields */
  FIELD_NPTISO_GV,
  /** index of \ref nptiso_struct::p_ext in \ref #fields */
  FIELD_NPTISO_PEXT,
  /** index of \ref nptiso_struct::p_inst in \ref #fields */
  FIELD_NPTISO_PINST,
  /** index of \ref nptiso_struct::p_inst_av in \ref #fields */
  FIELD_NPTISO_PINSTAV,
  /** index of \ref nptiso_struct::p_diff in \ref #fields */
  FIELD_NPTISO_PDIFF,
  /** index of \ref nptiso_struct::piston in \ref #fields */
  FIELD_NPTISO_PISTON,
  /** index of \ref #periodic in \ref #fields */
  FIELD_PERIODIC,
  /** index of \ref #skin in \ref #fields */
  FIELD_SKIN,
  /** index of \ref #temperature in \ref #fields */
  FIELD_TEMPERATURE,
  /** index of \ref thermo_switch in \ref #fields */
  FIELD_THERMO_SWITCH,
  /** index of \ref sim_time in  \ref #fields */
  FIELD_SIMTIME,
  /** index of \ref time_step in \ref #fields */
  FIELD_TIMESTEP,
  /** index of \ref timing_samples in  \ref #fields */
  FIELD_TIMINGSAMP,
  /** index of \ref max_cut_nonbonded in \ref #fields */
  FIELD_MCUT_NONBONDED,
  /** index of \ref verlet_reuse in  \ref #fields */
  FIELD_VERLETREUSE,
  /** index of \ref lattice_switch in \ref #fields */
  FIELD_LATTICE_SWITCH,
  /** index of \ref max_cut_bonded in \ref #fields */
  FIELD_MCUT_BONDED,
  /** index of \ref min_global_cut in \ref #fields */
  FIELD_MIN_GLOBAL_CUT,
  /** index of \ref ghmc_nmd in \ref #fields */
  FIELD_GHMC_NMD,
  /** index of \ref ghmc_phi in \ref #fields */
  FIELD_GHMC_PHI,
  /** index of \ref ghmc_phi in \ref #fields */
  FIELD_GHMC_RES,
  /** index of \ref ghmc_phi in \ref #fields */
  FIELD_GHMC_FLIP,
  /** index of \ref ghmc_phi in \ref #fields */
  FIELD_GHMC_SCALE,
  /** index of \ref lb_components in \ref #fields */
  FIELD_LB_COMPONENTS,
  /** index of \ref warnings in \ref #fields */
  FIELD_WARNINGS,
  /** index of \ref smaller_timestep in \ref #fields */
  FIELD_SMALLERTIMESTEP,
  /** index of \ref langevin_trans in \ref #fields */
  FIELD_LANGEVIN_TRANS_SWITCH,
  /** index of \ref langevin_rotate in \ref #fields */
  FIELD_LANGEVIN_ROT_SWITCH,
  /** index of \ref langevin_gamma_rotation in  \ref #fields */
  FIELD_LANGEVIN_GAMMA_ROTATION,
  FIELD_FORCE_CAP
};
/*@}*/

/** bool: whether to write out warnings or not */
extern int warnings;

#endif
