/*
  Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
  
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
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "virtual_sites_relative.hpp"
#include "collision.hpp"
#include "virtual_sites.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "communication.hpp" 
#include "parser.hpp" 

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
      Tcl_AppendResult(interp, s + 1, (char*) NULL);
    }

    if (collision_params.mode & COLLISION_MODE_VS) {
      sprintf(s, " bind_at_point_of_collision %f %d %d %d",
	      collision_params.distance, collision_params.bond_centers,
	      collision_params.bond_vs, collision_params.vs_particle_type);
      Tcl_AppendResult(interp, s + 1, (char*) NULL);
      return TCL_OK;
    }

    if (collision_params.mode & COLLISION_MODE_BIND_THREE_PARTICLES) {
      sprintf(s, " bind_three_particles %f %d %d %d",
	      collision_params.distance, collision_params.bond_centers,
	      collision_params.bond_three_particles, collision_params.three_particle_angle_resolution);
      Tcl_AppendResult(interp, s + 1, (char*) NULL);
    }
    else if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF) {
      sprintf(s, " glue_to_surface %f %d %d %d %d %d %d %f",
	      collision_params.distance, collision_params.bond_centers,
	      collision_params.bond_vs, collision_params.vs_particle_type,
	      collision_params.part_type_to_be_glued, 
	      collision_params.part_type_to_attach_vs_to,
	      collision_params.part_type_after_glueing,
	      collision_params.dist_glued_part_to_vs);
      Tcl_AppendResult(interp, s + 1, (char*) NULL);
    }
    else if (collision_params.mode & COLLISION_MODE_BOND) {
      sprintf(s, " bind_centers %f %d", collision_params.distance,
	      collision_params.bond_centers);
      Tcl_AppendResult(interp, s + 1, (char*) NULL);
    }
    // first character is always the separating space
    return TCL_OK;
  }

  argc--; argv++;

  // Otherwise, we set parameters
  if (ARG0_IS_S("off")) {
    collision_detection_set_params(0,0,0,0,0,0,0,0,0,0,0);
    return TCL_OK;
  }
  else {
    /* parameters of collision_detection_set_params */
    int mode = 0;
    
    // Distances
    double d,d2 = 0;

    // Bond types
    int bond_centers = 0;
    int bond_vs = 0;
    
    // Particle types for virtual sites based based methods
    int t,tg,tv,ta = 0;
    
    // /bond types for three particle binding
    int bond_three_particles=0;
    int angle_resolution=0;

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
    else if (ARG0_IS_S("glue_to_surface")) {
      mode |= COLLISION_MODE_BOND | COLLISION_MODE_GLUE_TO_SURF;
      if (argc != 9) {
	Tcl_AppendResult(interp, "Not enough parameters, need a distance, two bond types, four particle types and another distance as args.", (char*) NULL);
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
      if (!ARG_IS_I(5,tg)) {
	Tcl_AppendResult(interp, "Need a particle type as 5th arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(6,tv)) {
	Tcl_AppendResult(interp, "Need a particle type as 6th arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(7,ta)) {
	Tcl_AppendResult(interp, "Need a particle type as 7th arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_D(8,d2)) {
	Tcl_AppendResult(interp, "Need a distance as 8th arg.", (char*) NULL);
	return TCL_ERROR;
      }
      argc -= 9; argv += 8;
    }
    else if (ARG0_IS_S("bind_three_particles")) {
      mode |= COLLISION_MODE_BIND_THREE_PARTICLES | COLLISION_MODE_BOND;
      if (argc != 5) {
	Tcl_AppendResult(interp, "Not enough parameters, need a distance and two bond types.", (char*) NULL);
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
      if (!ARG_IS_I(3,bond_three_particles)) {
	Tcl_AppendResult(interp, "Need a bond type as 3rd arg.", (char*) NULL);
	return TCL_ERROR;
      }
      if (!ARG_IS_I(4,angle_resolution)) {
	Tcl_AppendResult(interp, "Need an angle resolution as 4th arg.", (char*) NULL);
	return TCL_ERROR;
      }
      argc -= 5; argv += 5;
    }
    else {
      Tcl_AppendResult(interp, "\"", argv[0], "\" is not a valid collision detection mode.", (char*) NULL);
      return TCL_ERROR;
    }
    
    int res = collision_detection_set_params(mode,d,bond_centers,bond_vs,t,d2,tg,tv,ta,bond_three_particles,angle_resolution);

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
    case 6:
      Tcl_AppendResult(interp, "Not enough angular bonds.", (char*) NULL);
      return TCL_ERROR;
    case 7:
      Tcl_AppendResult(interp, "bond_three_particles needs triple bonds.", (char*) NULL);
      return TCL_ERROR;
    }

    return TCL_OK;
  }
}

#endif
