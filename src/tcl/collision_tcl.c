/*
  Copyright (C) 2011,2012,2013 The ESPResSo project
  
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
#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "collision.h"
#include "virtual_sites.h"
#include "integrate.h"
#include "cells.h"
#include "communication.h" 
#include "parser.h" 

#ifdef COLLISION_DETECTION

int tclcommand_on_collision(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  // If no argumens are given, print status
  if (argc==1) {
    char s[128 + 3*TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE];

    if (collision_params.mode == 0) {
      sprintf(s, "off");
      Tcl_AppendResult(interp, s, (char*) NULL);
      return TCL_OK;
    }

    /* this one can be combined with the rest */
    if (collision_params.mode & COLLISION_MODE_EXCEPTION) {
      sprintf(s, " exception");
    }

    if (collision_params.mode & COLLISION_MODE_VS) {
      sprintf(s, " bind_at_point_of_collision %f %d %d %d",
	      collision_params.distance, collision_params.bond_centers,
	      collision_params.bond_vs, collision_params.vs_particle_type);
    }
    else if (collision_params.mode & COLLISION_MODE_BOND) {
      sprintf(s, " bind_centers %f %d", collision_params.distance,
	      collision_params.bond_centers);
    }
    // first character is always the separating space
    Tcl_AppendResult(interp, s + 1, (char*) NULL);
    return TCL_OK;
  }

  argc--; argv++;

  // Otherwise, we set parameters
  if (ARG0_IS_S("off")) {
    collision_detection_set_params(0,0,0,0,0);
    return TCL_OK;
  }
  else {
    /* parameters of collision_detection_set_params */
    int mode = 0;
    double d = 0;
    int bond_centers = 0;
    int bond_vs = 0;
    int t = 0;

    if (ARG0_IS_S("exception")) {
      mode = COLLISION_MODE_EXCEPTION;
      argc--; argv++;
    }
    if (argc == 0) {
      Tcl_AppendResult(interp, "throwing exception without creating any bond is not possible.", (char*) NULL);
      return TCL_ERROR;      
    }
    if (ARG0_IS_S("bind_centers")) {
      mode |= COLLISION_MODE_BOND;
      if (argc != 3) {
	Tcl_AppendResult(interp, "Not enough parameters, need a distance and a bond type as args.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_D(1,d)) {
	Tcl_AppendResult(interp, "Need a distance as 1st arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(2,bond_centers)) {
	Tcl_AppendResult(interp, "Need a bond type as 2nd argument.", (char*) NULL);
	return TCL_ERROR;
      }
      argc -= 3; argv += 3;
    }
    else if (ARG0_IS_S("bind_at_point_of_collision")) {
      mode |= COLLISION_MODE_BOND | COLLISION_MODE_VS;
      if (argc != 5) {
	Tcl_AppendResult(interp, "Not enough parameters, need a distance, two bond types, and a particle type as args.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_D(1,d)) {
	Tcl_AppendResult(interp, "Need a distance as 1st arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(2,bond_centers)) {
	Tcl_AppendResult(interp, "Need a bond type as 2nd arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(3,bond_vs)) {
	Tcl_AppendResult(interp, "Need a bond type as 3rd arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(4,t)) {
	Tcl_AppendResult(interp, "Need a particle type as 4th arg.", (char*) NULL);
	return TCL_ERROR;
      }
      argc -= 5; argv += 5;
    }
    else {
      Tcl_AppendResult(interp, "\"", argv[0], "\" is not a valid collision detection mode.", (char*) NULL);
      return TCL_ERROR;
    }
    
    int res = collision_detection_set_params(mode,d,bond_centers,bond_vs,t);

    switch (res) {
    case 1:
      Tcl_AppendResult(interp, "This mode requires the VIRTUAL_SITES_RELATIVE feature to be compiled in.", (char*) NULL);
      return TCL_ERROR;
    case 2:
      Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
      return TCL_ERROR;
    case 3:
      Tcl_AppendResult(interp, "Bond type does not exist.", (char*) NULL);
      return TCL_ERROR;
    case 4:
      Tcl_AppendResult(interp, "Real particles' bond has to be a pair bond.", (char*) NULL);
      return TCL_ERROR;
    case 5:
      Tcl_AppendResult(interp, "Virtual particles need a pair bond or triple bond.", (char*) NULL);
      return TCL_ERROR;
    }

    return TCL_OK;
  }
}

#endif
