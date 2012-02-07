/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
  for (int i=0; i < NUM_FEATURES; i++) {
    Tcl_AppendResult(interp, "{ ", (char *) NULL);
    Tcl_AppendResult(interp, FEATURES[i], (char *) NULL);
    Tcl_AppendResult(interp, " } ", (char *) NULL);
  }
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
