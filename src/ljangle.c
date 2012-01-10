/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

#include <math.h>
#include "ljangle.h"
#include "interaction_data.h"
#include "parser.h"


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
    @param ljangleforcecap the maximal force, 0 to disable, -1 for individual cutoff
    for each of the interactions.
*/
int ljangleforcecap_set_params(double ljangleforcecap)
{
  if (ljangle_force_cap != -1.0)
    mpi_ljangle_cap_forces(ljangle_force_cap);
   
  return TCL_OK;
}

int ljangle_set_params(int part_type_a, int part_type_b,
				double eps, double sig, double cut,
				int b1p, int b1n, int b2p, int b2n,
				double cap_radius, double z0, double dz, 
				double kappa, double epsprime)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJ_ANGLE should be symmetrically */
  data->LJANGLE_eps         = data_sym->LJANGLE_eps         = eps;
  data->LJANGLE_sig         = data_sym->LJANGLE_sig         = sig;
  data->LJANGLE_cut         = data_sym->LJANGLE_cut         = cut;
  data->LJANGLE_bonded1pos  = data_sym->LJANGLE_bonded1pos  = b1p;
  data->LJANGLE_bonded1neg  = data_sym->LJANGLE_bonded1neg  = b1n;
  data->LJANGLE_bonded2pos  = data_sym->LJANGLE_bonded2pos  = b2p;
  data->LJANGLE_bonded2neg  = data_sym->LJANGLE_bonded2neg  = b2n;
  
  data->LJANGLE_bonded1type     =   data_sym->LJANGLE_bonded1type = part_type_a;

  if (cap_radius > 0) {
    data->LJANGLE_capradius = cap_radius;
    data_sym->LJANGLE_capradius = cap_radius;
  }

  if (dz > 0.) {
    data->LJANGLE_z0       = data_sym->LJANGLE_z0       = z0;
    data->LJANGLE_dz       = data_sym->LJANGLE_dz       = dz;
    data->LJANGLE_kappa    = data_sym->LJANGLE_kappa    = kappa;
    data->LJANGLE_epsprime = data_sym->LJANGLE_epsprime = epsprime;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (ljangle_force_cap != -1.0)
    mpi_ljangle_cap_forces(ljangle_force_cap);
    
  return TCL_OK;
}

/// parser for the forcecap
int tclcommand_inter_parse_ljangleforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE];
  if (argc == 0) {
    if (ljangle_force_cap == -1.0)
      Tcl_AppendResult(interp, "ljangleforcecap individual", (char *) NULL);
    else {
      Tcl_PrintDouble(interp, ljangle_force_cap, buffer);
      Tcl_AppendResult(interp, "ljangleforcecap ", buffer, (char *) NULL); 
    }
    return TCL_OK;
  }
  if (argc > 1) {
    Tcl_AppendResult(interp, "inter ljangleforcecap takes at most 1 parameter",
		     (char *) NULL);
    return TCL_ERROR;
  }
  if (ARG0_IS_S("individual"))
    ljangle_force_cap = -1.0;
  else if (! ARG0_IS_D(ljangle_force_cap) || ljangle_force_cap < 0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		     (char *) NULL);
    return TCL_ERROR;
  }
  CHECK_VALUE(ljangleforcecap_set_params(ljangle_force_cap),
	      "If you can read this, you should change it. (Use the source Luke!)");
  return TCL_ERROR; 
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


/** calculate ljangle_capradius from ljangle_force_cap */
/* This routine does not take the optional 2nd environment into account. */
void calc_ljangle_cap_radii(double force_cap)
{
  int i,j,cnt=0;
  IA_parameters *params;
  double force=0.0, rad=0.0, step, frac2, frac10;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      if(force_cap > 0.0 && params->LJANGLE_eps > 0.0) {
	/* I think we have to solve this numerically... and very crude as well */
	cnt=0;
	rad = params->LJANGLE_sig;
	step = -0.1 * params->LJANGLE_sig;
	force=0.0;
	
	while(step != 0) {
	  frac2  = SQR(params->LJANGLE_sig/rad);
	  frac10 = frac2*frac2*frac2*frac2*frac2;
	  force = 60.0 * params->LJANGLE_eps * frac10*(frac2 - 1.0) / rad;
	  if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
	    step = - (step/2.0); 
	  }
	  if(fabs(force-force_cap) < 1.0e-6) step=0;
	  rad += step; cnt++;
	} 
	params->LJANGLE_capradius = rad;
      }
      else {
	params->LJANGLE_capradius = 0.0; 
      }
      FORCE_TRACE(fprintf(stderr,"%d: LJANGLE Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
			  this_node,i,j,rad,force,cnt));
    }
  }
}



