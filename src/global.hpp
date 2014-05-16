/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

/** Field is of type integer in \ref Datafield. */
#define TYPE_INT    0
/** Field is of type double in \ref Datafield. */
#define TYPE_DOUBLE 1
/** Field is of type bool, i.e. bit array, in \ref Datafield.
    Note that the field is stored in whatever an integer is.
    I guess you can at least assume 16 bits...
*/
#define TYPE_BOOL 2

/** Maximal size of an array in \ref Datafield. */
#define MAX_DIMENSION 64

/** Type describing global variables. These are accessible from the
    front end, and are distributed to all compute nodes. */
typedef struct {
  /** Physical address of the variable. */
  void        *data;
  /** Type of the variable, either \ref TYPE_INT or \ref TYPE_DOUBLE.*/
  int          type;
  /** Dimension of the variable. Limited to \ref MAX_DIMENSION */
  int          dimension;
  /** Name of the variable, mainly used for the front end and debugging */
  const char  *name;
  /** Minimal number of characters needed for unique identification of the
      variable. */
  int min_length;
} Datafield;

/** This array contains the description of all global variables that
    are synchronized across nodes and that can be changed/adressed via
    the TCL command setmd. read the documentation of \ref Datafield
    befor you add new features. */
extern const Datafield fields[];

/** \name Field Enumeration
    These numbers identify the variables given in
    \ref #fields for use with \ref mpi_bcast_parameter.
*/
/*@{*/
/** index of \ref box_l in \ref #fields */
#define FIELD_BOXL                0  
/** index of \ref DomainDecomposition::cell_grid in  \ref #fields */
#define FIELD_CELLGRID            1
/** index of \ref DomainDecomposition::cell_size in  \ref #fields */
#define FIELD_CELLSIZE            2
/** index of \ref dpd_gamma in  \ref #fields */
#define FIELD_DPD_GAMMA           3
/** index of \ref dpd_r_cut in  \ref #fields */
#define FIELD_DPD_RCUT            4
/** index of \ref langevin_gamma in  \ref #fields */
#define FIELD_LANGEVIN_GAMMA      5
/** index of \ref integ_switch in \ref #fields */
#define FIELD_INTEG_SWITCH        6
/** index of \ref local_box_l in \ref #fields */
#define FIELD_LBOXL               7
/** index of \ref max_cut in \ref #fields */
#define FIELD_MCUT                8
/** index of \ref max_num_cells  in \ref #fields */
#define FIELD_MAXNUMCELLS         9
/** index of \ref max_seen_particle in \ref #fields */
#define FIELD_MAXPART             10
/** index of \ref max_range in \ref #fields */
#define FIELD_MAXRANGE            11
/** index of \ref max_skin in  \ref #fields */
#define FIELD_MAXSKIN             12
/** index of \ref min_num_cells  in \ref #fields */
#define FIELD_MINNUMCELLS         13
/** index of \ref n_layers in  \ref #fields */
#define FIELD_NLAYERS             14
/** index of \ref n_nodes in \ref #fields */
#define FIELD_NNODES              15
/** index of \ref n_part in  \ref #fields */
#define FIELD_NPART               16
/** index of \ref n_particle_types in \ref #fields */
#define FIELD_NPARTTYPE           17
/** index of \ref n_rigidbonds in \ref #fields */
#define FIELD_RIGIDBONDS          18
/** index of \ref node_grid in \ref #fields */
#define FIELD_NODEGRID            19
/** index of \ref nptiso_gamma0 in \ref #fields */
#define FIELD_NPTISO_G0           20
/** index of \ref nptiso_gammav in \ref #fields */
#define FIELD_NPTISO_GV           21
/** index of \ref nptiso_struct::p_ext in \ref #fields */
#define FIELD_NPTISO_PEXT         22      
/** index of \ref nptiso_struct::p_inst in \ref #fields */
#define FIELD_NPTISO_PINST        23
/** index of \ref nptiso_struct::p_inst_av in \ref #fields */
#define FIELD_NPTISO_PINSTAV      24
/** index of \ref nptiso_struct::p_diff in \ref #fields */
#define FIELD_NPTISO_PDIFF        25
/** index of \ref nptiso_struct::piston in \ref #fields */
#define FIELD_NPTISO_PISTON       26
/** index of \ref #periodic in \ref #fields */
#define FIELD_PERIODIC            27
/** index of \ref #skin in \ref #fields */
#define FIELD_SKIN                28
/** index of \ref #temperature in \ref #fields */
#define FIELD_TEMPERATURE         29
/** index of \ref thermo_switch in \ref #fields */
#define FIELD_THERMO_SWITCH       30
/** index of \ref sim_time in  \ref #fields */
#define FIELD_SIMTIME             31
/** index of \ref time_step in \ref #fields */
#define FIELD_TIMESTEP            32
/** index of \ref timing_samples in  \ref #fields */
#define FIELD_TIMINGSAMP          33
/** index of \ref transfer_rate  in \ref #fields */
#define FIELD_TRANSFERRATE        34
/** index of \ref max_cut_nonbonded in \ref #fields */
#define FIELD_MCUT_NONBONDED      35
/** index of \ref verlet_reuse in  \ref #fields */
#define FIELD_VERLETREUSE         36
/** index of \ref lattice_switch in \ref #fields */
#define FIELD_LATTICE_SWITCH      37
/** index of \ref dpd_tgamma in \ref #fields */
#define FIELD_DPD_TGAMMA          38
/** index of \ref dpd_tr_cut in \ref #fields */
#define FIELD_DPD_TRCUT           39
/** index of \ref dpd_twf in \ref #fields */
#define FIELD_DPD_TWF             40
/** index of \ref dpd_wf in \ref #fields */
#define FIELD_DPD_WF              41
/** index of \ref max_cut_bonded in \ref #fields */
#define FIELD_MCUT_BONDED         42
/** index of \ref transfer_rate in \ref #fields */
#define FIELD_TRANSFER_RATE       43
/** index of \ref min_global_cut in \ref #fields */
#define FIELD_MIN_GLOBAL_CUT      44
/** index of \ref ghmc_nmd in \ref #fields */
#define FIELD_GHMC_NMD            45
/** index of \ref ghmc_phi in \ref #fields */
#define FIELD_GHMC_PHI            46
/** index of \ref ghmc_phi in \ref #fields */
#define FIELD_GHMC_RES            47 
/** index of \ref ghmc_phi in \ref #fields */
#define FIELD_GHMC_FLIP           48
/** index of \ref ghmc_phi in \ref #fields */
#define FIELD_GHMC_SCALE          49 
/** index of \ref lb_components in \ref #fields */
#define FIELD_LB_COMPONENTS       50 
/** index of \ref warnings in \ref #fields */
#define FIELD_WARNINGS            51
/** DPD_IGNORE_FIXED_PARTICLES */
#define FIELD_DPD_IGNORE_FIXED_PARTICLES 52
/*@}*/

/** bool: whether to write out warnings or not */
extern int warnings;

#endif
