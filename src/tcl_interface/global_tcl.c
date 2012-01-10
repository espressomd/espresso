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
/** \file global.c
    Implementation of \ref global.h "global.h".
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "global_tcl.h"
/* from these modules we modify variables: */
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "thermostat.h"
#include "forces.h"
#include "verlet.h"
#include "p3m.h"
#include "imd.h"
#include "tuning.h"
#include "domain_decomposition.h"
#include "layered.h"
#include "pressure.h"
#include "rattle.h"
#include "lattice.h"
#include "adresso.h"
#include "tcl_interface/integrate_tcl.h"

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

int tclcommand_code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  if (argc < 2) {
    tclcallback_version(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    tclcallback_compilation(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    tclcallback_debug(interp);
  }
  else {
    if(!strncmp(argv[1], "version" , strlen(argv[1]))) {
      tclcallback_version(interp);
    }
    else if(!strncmp(argv[1], "compilation" , strlen(argv[1]))) {
      tclcallback_compilation(interp);
    }
    else if(!strncmp(argv[1], "debug" , strlen(argv[1]))) {
      tclcallback_debug(interp);
    }
    else {
      Tcl_AppendResult(interp, "info ",argv[1]," not known!", (char *) NULL);
    }
  }
  return (TCL_OK);
}
