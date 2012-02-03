/*
  Copyright (C) 2010,2011 The ESPResSo project
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
/** \file config_tcl.c
 *
 *  contains code_info and version stuff.
*/
#include "utils.h"
#include <tcl.h>

static int tclcommand_code_info_version(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, PACKAGE_NAME ": " PACKAGE_VERSION, (char *) NULL);
  return (TCL_OK);
}

/** callback for compilation status. */
static int tclcommand_code_info_compilation(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Compilation status ", (char *) NULL);
#ifdef FFTW
  Tcl_AppendResult(interp, "{ FFTW } ", (char *) NULL);
#endif
#ifdef CUDA
  Tcl_AppendResult(interp, "{ CUDA } ", (char *) NULL);
#endif
#ifdef TK
  Tcl_AppendResult(interp, "{ TK } ", (char *) NULL);
#endif
#ifdef MODES
  Tcl_AppendResult(interp, "{ MODES } ", (char *) NULL);
#endif
#ifdef PARTIAL_PERIODIC
  Tcl_AppendResult(interp, "{ PARTIAL_PERIODIC } ", (char *) NULL);
#endif
#ifdef LJ_WARN_WHEN_CLOSE
  Tcl_AppendResult(interp, "{ LJ_WARN_WHEN_CLOSE } ", (char *) NULL);
#endif
#ifdef ELECTROSTATICS
  Tcl_AppendResult(interp, "{ ELECTROSTATICS } ", (char *) NULL);
#endif
#ifdef ROTATION
  Tcl_AppendResult(interp, "{ ROTATION } ", (char *) NULL);
#endif
#ifdef ROTATIONAL_INERTIA
  Tcl_AppendResult(interp, "{ ROTATIONAL_INERTIA } ", (char *) NULL);
#endif
#ifdef DIPOLES
  Tcl_AppendResult(interp, "{ DIPOLES } ", (char *) NULL);
#endif
#ifdef MASS
  Tcl_AppendResult(interp, "{ MASS } ", (char *) NULL);
#endif
#ifdef EXTERNAL_FORCES
  Tcl_AppendResult(interp, "{ EXTERNAL_FORCES } ", (char *) NULL);
#endif
#ifdef CONSTRAINTS
  Tcl_AppendResult(interp, "{ CONSTRAINTS } ", (char *) NULL);
#endif
#ifdef BOND_CONSTRAINT
  Tcl_AppendResult(interp, "{ BOND_CONSTRAINT } ", (char *) NULL);
#endif
#ifdef BOND_VIRTUAL
  Tcl_AppendResult(interp, "{ BOND_VIRTUAL } ", (char *) NULL);
#endif
#ifdef COLLISION_DETECTION 
  Tcl_AppendResult(interp, "{ COLLISION_DETECTION } ", (char *) NULL);
#endif
#ifdef EXCLUSIONS
  Tcl_AppendResult(interp, "{ EXCLUSIONS } ", (char *) NULL);
#endif
#ifdef COMFORCE
  Tcl_AppendResult(interp, "{ COMFORCE } ", (char *) NULL);
#endif
#ifdef COMFIXED
  Tcl_AppendResult(interp, "{ COMFIXED } ", (char *) NULL);
#endif
#ifdef TABULATED
  Tcl_AppendResult(interp, "{ TABULATED } ", (char *) NULL);
#endif
#ifdef OVERLAPPED 
  Tcl_AppendResult(interp, "{ OVERLAPPED } ", (char *) NULL);
#endif
#ifdef LENNARD_JONES
  Tcl_AppendResult(interp, "{ LENNARD_JONES } ", (char *) NULL);
#endif
#ifdef LENNARD_JONES_GENERIC
  Tcl_AppendResult(interp, "{ LENNARD_JONES_GENERIC } ", (char *) NULL);
#endif
#ifdef GAY_BERNE
  Tcl_AppendResult(interp, "{ GAY_BERNE } ", (char *) NULL);
#endif
#ifdef LJ_ANGLE
  Tcl_AppendResult(interp, "{ LJ_ANGLE } ", (char *) NULL);
#endif
#ifdef SMOOTH_STEP
  Tcl_AppendResult(interp, "{ SMOOTH_STEP } ", (char *) NULL);
#endif
#ifdef HERTZIAN
  Tcl_AppendResult(interp, "{ HERTZIAN } ", (char *) NULL);
#endif
#ifdef BMHTF_NACL
  Tcl_AppendResult(interp, "{ BMHTF_NACL } ", (char *) NULL);
#endif
#ifdef MORSE
  Tcl_AppendResult(interp, "{ MORSE } ", (char *) NULL);
#endif
#ifdef BUCKINGHAM
  Tcl_AppendResult(interp, "{ BUCKINGHAM } ", (char *) NULL);
#endif
#ifdef SOFT_SPHERE
  Tcl_AppendResult(interp, "{ SOFT_SPHERE } ", (char *) NULL);
#endif
#ifdef LJCOS
  Tcl_AppendResult(interp, "{ LJCOS } ", (char *) NULL);
#endif
#ifdef LJCOS2
  Tcl_AppendResult(interp, "{ LJCOS2 } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_HARMONIC
  Tcl_AppendResult(interp, "{ BOND_ANGLE_HARMONIC } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_COSINE
  Tcl_AppendResult(interp, "{ BOND_ANGLE_COSINE } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_COSSQUARE
  Tcl_AppendResult(interp, "{ BOND_ANGLE_COSSQUARE } ", (char *) NULL);
#endif
#ifdef BOND_ANGLEDIST_HARMONIC
  Tcl_AppendResult(interp, "{ BOND_ANGLEDIST_HARMONIC } ", (char *) NULL);
#endif
#ifdef BOND_ENDANGLEDIST_HARMONIC
  Tcl_AppendResult(interp, "{ BOND_ENDANGLEDIST_HARMONIC } ", (char *) NULL);
#endif
#ifdef NEMD
  Tcl_AppendResult(interp, "{ NEMD } ", (char *) NULL);
#endif
#ifdef NPT
  Tcl_AppendResult(interp, "{ NPT } ", (char *) NULL);
#endif
#ifdef LB_GPU
  Tcl_AppendResult(interp, "{ LB_GPU } ", (char *) NULL);
#endif
#ifdef TRANS_DPD
  Tcl_AppendResult(interp, "{ TRANS_DPD } ", (char *) NULL);
#endif
#ifdef DPD
  Tcl_AppendResult(interp, "{ DPD } ", (char *) NULL);
#endif
#ifdef DPD_MASS_RED
  Tcl_AppendResult(interp, "{ DPD_MASS_RED } ", (char *) NULL);
#endif
#ifdef DPD_MASS_LIN
  Tcl_AppendResult(interp, "{ DPD_MASS_LIN } ", (char *) NULL);
#endif
#ifdef LB
  Tcl_AppendResult(interp, "{ LB } ", (char *) NULL);
#endif
#ifdef LB_BOUNDARIES
  Tcl_AppendResult(interp, "{ LB_BOUNDARIES } ", (char *) NULL);
#endif
#ifdef LB_BOUNDARIES_GPU
  Tcl_AppendResult(interp, "{ LB_BOUNDARIES_GPU } ", (char *) NULL);
#endif
#ifdef INTER_DPD
  Tcl_AppendResult(interp, "{ INTER_DPD } ", (char *) NULL);
#endif
#ifdef INTER_RF
  Tcl_AppendResult(interp, "{ INTER_RF } ", (char *) NULL);
#endif
#ifdef NO_INTRA_NB
  Tcl_AppendResult(interp, "{ NO_INTRA_NB } ", (char *) NULL);
#endif
#ifdef VIRTUAL_SITES_COM
  Tcl_AppendResult(interp, "{ VIRTUAL_SITES_COM } ", (char *) NULL);
#endif
#ifdef VIRTUAL_SITES_RELATIVE
  Tcl_AppendResult(interp, "{ VIRTUAL_SITES_RELATIVE } ", (char *) NULL);
#endif
#ifdef THERMOSTAT_IGNORE_NON_VIRTUAL
  Tcl_AppendResult(interp, "{ THERMOSTAT_IGNORE_NON_VIRTUAL } ", (char *) NULL);
#endif
#ifdef METADYNAMICS
  Tcl_AppendResult(interp, "{ METADYNAMICS } ", (char *) NULL);
#endif
#ifdef MOL_CUT
  Tcl_AppendResult(interp, "{ MOL_CUT } ", (char *) NULL);
#endif
#ifdef TUNABLE_SLIP
  Tcl_AppendResult(interp, "{ TUNABLE_SLIP } ", (char *) NULL);
#endif
#ifdef ADRESS
  Tcl_AppendResult(interp, "{ ADRESS } ", (char *) NULL);
#ifdef ADRESS_INIT
  Tcl_AppendResult(interp, "{ ADRESS_INIT } ", (char *) NULL);
#endif
#ifdef INTERFACE_CORRECTION
  Tcl_AppendResult(interp, "{ INTERFACE_CORRECTION } ", (char *) NULL);
#endif
#ifdef ROTATION_PER_PARTICLE
  Tcl_AppendResult(interp, "{ ROTATION_PER_PARTICLE } ", (char *) NULL);
#endif
#ifdef LANGEVIN_PER_PARTICLE
  Tcl_AppendResult(interp, "{ LANGEVIN_PER_PARTICLE } ", (char *) NULL);
#endif
#endif
 Tcl_AppendResult(interp, "}", (char *) NULL);
  return (TCL_OK);
}

static int tclcommand_code_info_debug(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Debug status { ", (char *) NULL);
#ifdef COMM_DEBUG
  Tcl_AppendResult(interp, "COMM_DEBUG ", (char *) NULL);
#endif
#ifdef INTEG_DEBUG
  Tcl_AppendResult(interp, "INTEG_DEBUG ", (char *) NULL);
#endif
#ifdef CELL_DEBUG
  Tcl_AppendResult(interp, "CELL_DEBUG ", (char *) NULL);
#endif
#ifdef GHOST_DEBUG
  Tcl_AppendResult(interp, "GHOST_DEBUG ", (char *) NULL);
#endif
#ifdef GRID_DEBUG 
  Tcl_AppendResult(interp, "GRID_DEBUG ", (char *) NULL);
#endif
#ifdef VERLET_DEBUG
  Tcl_AppendResult(interp, "VERLET_DEBUG ", (char *) NULL);
#endif
#ifdef PARTICLE_DEBUG
  Tcl_AppendResult(interp, "PARTICLE_DEBUG ", (char *) NULL);
#endif
#ifdef P3M_DEBUG
  Tcl_AppendResult(interp, "P3M_DEBUG ", (char *) NULL);
#endif
#ifdef FFT_DEBUG
  Tcl_AppendResult(interp, "FFT_DEBUG ", (char *) NULL);
#endif
#ifdef RANDOM_DEBUG
  Tcl_AppendResult(interp, "RANDOM_DEBUG ", (char *) NULL);
#endif
#ifdef FORCE_DEBUG
  Tcl_AppendResult(interp, "FORCE_DEBUG ", (char *) NULL);
#endif
#ifdef THERMO_DEBUG
  Tcl_AppendResult(interp, "THERMO_DEBUG ", (char *) NULL);
#endif
#ifdef LJ_DEBUG
  Tcl_AppendResult(interp, "LJ_DEBUG ", (char *) NULL);
#endif
#ifdef ESR_DEBUG
  Tcl_AppendResult(interp, "ESR_DEBUG ", (char *) NULL);
#endif
#ifdef ESK_DEBUG
  Tcl_AppendResult(interp, "ESK_DEBUG ", (char *) NULL);
#endif
#ifdef FENE_DEBUG
  Tcl_AppendResult(interp, "FENE_DEBUG ", (char *) NULL);
#endif
#ifdef GHOST_FORCE_DEBUG
  Tcl_AppendResult(interp, "GHOST_FORCE_DEBUG ", (char *) NULL);
#endif
#ifdef ONEPART_DEBUG
  Tcl_AppendResult(interp, "ONEPART_DEBUG ", (char *) NULL);
#endif
#ifdef STAT_DEBUG
  Tcl_AppendResult(interp, "STAT_DEBUG ", (char *) NULL);
#endif
#ifdef POLY_DEBUG
  Tcl_AppendResult(interp, "POLY_DEBUG ", (char *) NULL);
#endif
#ifdef MPI_CORE
  Tcl_AppendResult(interp, "MPI_CORE ", (char *) NULL);
#endif
#ifdef FORCE_CORE
  Tcl_AppendResult(interp, "FORCE_CORE ", (char *) NULL);
#endif
#ifdef ADDITIONAL_CHECKS
  Tcl_AppendResult(interp, "ADDITIONAL_CHECKS ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "} }", (char *) NULL);
  return (TCL_OK);
}

int tclcommand_code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  if (argc < 2) {
    tclcommand_code_info_version(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    tclcommand_code_info_compilation(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    tclcommand_code_info_debug(interp);
  }
  else {
    if(!strncmp(argv[1], "version" , strlen(argv[1]))) {
      tclcommand_code_info_version(interp);
    }
    else if(!strncmp(argv[1], "compilation" , strlen(argv[1]))) {
      tclcommand_code_info_compilation(interp);
    }
    else if(!strncmp(argv[1], "debug" , strlen(argv[1]))) {
      tclcommand_code_info_debug(interp);
    }
    else {
      Tcl_AppendResult(interp, "info ",argv[1]," not known!", (char *) NULL);
    }
  }
  return (TCL_OK);
}
