/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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

/** \file ljangle.h
 *  Routines to calculate the lennard-jones 12-10 with angular dependance.
 *  The potential is a product of a 12-10 LJ potential with two cos^2.
 *  The potential actually relies on 6 particles: the 2 primary beads
 *  and each bead needs two other particles to define an orientation.
 *  We calculate the interaction explicitly *without* the use of the ROTATION feature.
 *
 *  Optional: simulate two different environments in which the interaction
 *  strengths differ. For example: simulate hydrogen-bonds both in water and
 *  inside a bilayer. The two environments are distinguished by their
 *  z-position. Input: the midplane of the 2nd environment, its total
 *  thickness, the thickness of the interface, and the value of the
 *  interaction strength in this medium. The interaction strengh of the second
 *  environment must be *stronger* than of the first one.
 *
 *  \ref forces.c
 */

#include "utils.hpp"

#ifdef LJ_ANGLE
#include "ljangle.hpp"
#include "ljangle_tcl.hpp"
#include <cmath>
#include "interaction_data.hpp"
#include "parser.hpp"
#include "communication.hpp"
#include "forcecap_tcl.hpp"

int tclprint_to_result_ljangleIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->LJANGLE_eps, buffer);
  Tcl_AppendResult(interp, "LJ-Angle ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  // Who are the bonded partners?
  Tcl_PrintDouble(interp, data->LJANGLE_bonded1pos, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_bonded1neg, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_bonded2pos, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_bonded2neg, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  // Optional argument: cap radius
  Tcl_PrintDouble(interp, data->LJANGLE_capradius, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  // Optional arguments: simulate two different interaction strengths
  Tcl_PrintDouble(interp, data->LJANGLE_z0, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_dz, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_kappa, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJANGLE_epsprime, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  
  return TCL_OK;
}

/** set the force cap for the directional LJ interaction.
    @param interp tcl interpreter handle
    @param argc number of arguments.
    @param argv arguments.
*/

int tclcommand_inter_parse_ljangleforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  if (argc==1) {
    fprintf(stderr, "WARNING: \"inter ljangleforcecap\" is deprecated "
	    "and will be removed in some further version. "
	    "Use \"inter forcecap\" instead.\n");
  }
  return tclcommand_inter_parse_forcecap(interp, argc, argv);
}



int tclcommand_inter_parse_ljangle(Tcl_Interp * interp,
			    int part_type_a, int part_type_b,
			    int argc, char ** argv)
{
  /* parameters needed for LJANGLE */
  double eps, sig, cut, cap_radius, z0, dz, kappa, epsprime;
  int b1p, b1n, b2p, b2n;
  int change;


  /* get lj-hbond interaction type */
  if (argc < 8) {
    Tcl_AppendResult(interp, "lj-angle needs 7 parameters: "
		     "<ljangle_eps> <ljangle_sig> <ljangle_cut> "
		     "<ljangle_1st_bonded_pos> <ljangle_1st_bonded_neg> "
		     "<ljangle_2nd_bonded_pos> <ljangle_2nd_bonded_neg>",
		     (char *) NULL);
    return 0;
  }

  /* copy lj-angle parameters */
  if ((! ARG_IS_D(1, eps))    ||
      (! ARG_IS_D(2, sig))    ||
      (! ARG_IS_D(3, cut))    ||
      (! ARG_IS_I(4, b1p))    ||
      (! ARG_IS_I(5, b1n))    ||
      (! ARG_IS_I(6, b2p))    ||
      (! ARG_IS_I(7, b2n))) {
    Tcl_AppendResult(interp, "lj-angle needs 3 DOUBLE  and 4 INT parameters: "
		     "<ljangle_eps> <ljangle_sig> <ljangle_cut> "
		     "<ljangle_1st_bonded_pos> <ljangle_1st_bonded_neg> "
		     "<ljangle_2nd_bonded_pos> <ljangle_2nd_bonded_neg>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (b1p==b1n || b2p==b2n) {
    Tcl_AppendResult(interp,"lj-angle needs 2 *different* bonded partners for each particle.",(char *) NULL);
    return TCL_ERROR;
  }

  if (b1p==0 || b1n==0 || b2p==0 || b2n==0) {
    Tcl_AppendResult(interp,"lj-angle: one of the bonded partners *is* the particle itself.",(char *) NULL);
    return TCL_ERROR;
  }  

  change = 8;
	
  cap_radius = -1.0;
  z0 = 0.;
  dz = -1;
  kappa = 0.;
  epsprime = 0.;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 9 && ARG_IS_D(8, cap_radius)) {
    change++;
    /* Check for the optional parameters describing the second environment */
    if (argc == 13 && ARG_IS_D(9, z0) && ARG_IS_D(10, dz) && ARG_IS_D(11, kappa) && ARG_IS_D(12, epsprime)) {
      /* check the integrity of the variables */
      if (z0 < 0. || z0 > box_l[2] || dz < 0. || dz > box_l[2] || kappa < 0. || kappa >= dz) {
	Tcl_AppendResult(interp,"lj-angle: Optional parameters for 2nd environment are not compatible with the box size.",(char *) NULL);
	return TCL_ERROR;
      }
      change += 4;
    }
    else
      Tcl_ResetResult(interp);
  }
  else
    Tcl_ResetResult(interp);
  if (ljangle_set_params(part_type_a, part_type_b,
			 eps, sig, cut, b1p, b1n, b2p, b2n, 
			 cap_radius, z0, dz, kappa, epsprime) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  
  return change;
}

#endif
