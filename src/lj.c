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

#include "lj.h"
#include "parser.h"


#ifdef LENNARD_JONES
#include "mol_cut.h"

int tclprint_to_result_ljIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->LJ_eps, buffer);
  Tcl_AppendResult(interp, "lennard-jones ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJ_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJ_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJ_shift, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJ_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJ_capradius, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  Tcl_PrintDouble(interp, data->LJ_min, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  

  
  return TCL_OK;
}

/** set the force cap for the LJ interaction.
    @param ljforcecap the maximal force, 0 to disable, -1 for individual cutoff
    for each of the interactions.
*/
int ljforcecap_set_params(double ljforcecap)
{
  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);
  
  return TCL_OK;
}

int lennard_jones_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double cut,
				      double shift, double offset,
				      double cap_radius, double min)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJ should be symmetrically */
  data->LJ_eps    = data_sym->LJ_eps    = eps;
  data->LJ_sig    = data_sym->LJ_sig    = sig;
  data->LJ_cut    = data_sym->LJ_cut    = cut;
  data->LJ_shift  = data_sym->LJ_shift  = shift;
  data->LJ_offset = data_sym->LJ_offset = offset;
 
  if (cap_radius > 0) {
    data->LJ_capradius = cap_radius;
    data_sym->LJ_capradius = cap_radius;
  }

  if (min > 0) {
	  data->LJ_min = min;
	  data_sym->LJ_min = min;	  
  }
  

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}

/// parser for the forcecap
int tclcommand_inter_parse_ljforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE];

  if (argc == 0) {
    if (lj_force_cap == -1.0)
      Tcl_AppendResult(interp, "ljforcecap individual", (char *) NULL);
    else {
      Tcl_PrintDouble(interp, lj_force_cap, buffer);
      Tcl_AppendResult(interp, "ljforcecap ", buffer, (char *) NULL);
    }
    return TCL_OK;
  }

  if (argc > 1) {
    Tcl_AppendResult(interp, "inter ljforcecap takes at most 1 parameter",
		     (char *) NULL);      
    return TCL_ERROR;
  }
  
  if (ARG0_IS_S("individual"))
      lj_force_cap = -1.0;
  else if (! ARG0_IS_D(lj_force_cap) || lj_force_cap < 0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(ljforcecap_set_params(lj_force_cap),
	      "If you can read this, you should change it. (Use the source Luke!)");
  return TCL_ERROR;
}

int tclcommand_inter_parse_lj(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJ */
  double eps, sig, cut, shift, offset, cap_radius, min;
  int compute_auto_shift, change;

  /* get lennard-jones interaction type */
  if (argc < 4) {
    Tcl_AppendResult(interp, "lennard-jones needs at least 3 parameters: "
		     "<lj_eps> <lj_sig> <lj_cut> [<lj_shift> [<lj_offset> [<lj_cap> [<lj_min>]]]]",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 1;

  /* PARSE LENNARD-JONES COMMAND LINE */
  /* epsilon */
  if (! ARG_IS_D(1, eps)) {
    Tcl_AppendResult(interp, "<lj_eps> should be a double",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change++;

  /* sigma */
  if (! ARG_IS_D(2, sig)) {
    Tcl_AppendResult(interp, "<lj_sig> should be a double",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change++;

  /* cutoff */
  if (! ARG_IS_D(3, cut)) {
    Tcl_AppendResult(interp, "<lj_cut> should be a double",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change++;
  
  /* shift */
  if (argc > 4) {
    if (ARG_IS_D(4, shift)) {
      /* if a shift is given, use that one */
      compute_auto_shift = 0;
      change++;
    } else if (ARG_IS_S(4, "auto")) {
      compute_auto_shift = 1;
      /* if shift is "auto", autocompute the shift */
      change++;
    } else {
      Tcl_AppendResult(interp, "<lj_shift> should be a double or \"auto\"",
		       (char *) NULL);
      return TCL_ERROR;
    }
  } else {
    /* by default, compute the shift automatically */
    compute_auto_shift = 1;
  }
  /* the shift can be computed automatically only when the other
     parameters have been determined, see below */

  /* offset */
  if (argc > 5) {
    if (!ARG_IS_D(5, offset)) {
      Tcl_AppendResult(interp, "<lj_off> should be a double",
		       (char *) NULL);
      return TCL_ERROR;
    }
    change++;
  } else {
    offset = 0.0;
  }
  
  /* cap_radius */
  if (argc > 6) {
    if (!ARG_IS_D(6, cap_radius)) {
      Tcl_AppendResult(interp, "<lj_cap> should be a double",
		       (char *) NULL);
      return TCL_ERROR;
    }
    change++;
  } else {
    cap_radius = -1.0;
  }

  /* minimal radius */
  if (argc > 7) {
    if (!ARG_IS_D(7, min)) {
      Tcl_AppendResult(interp, "<lj_rmin> should be a double",
		       (char *) NULL);
      return TCL_ERROR;
    }
    change ++;
  } else {
    min = 0.0;
  }

  /* automatically compute the shift */
  if (compute_auto_shift)
    shift = -(pow(sig/cut, 12) - pow(sig/cut, 6));

  /* now set the parameters */
  if (lennard_jones_set_params(part_type_a, part_type_b,
			       eps, sig, cut, shift, offset,
			       cap_radius, min) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", 
		     (char *) NULL);
    return 0;
  }
  return change;
}


/** Calculate lennard Jones force between particle p1 and p2 */


/** calculate lj_capradius from lj_force_cap */
void calc_lj_cap_radii(double force_cap)
{
  int i,j,cnt=0;
  IA_parameters *params;
  double force=0.0, rad=0.0, step, frac2, frac6;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      if(force_cap > 0.0 && params->LJ_eps > 0.0) {
	/* I think we have to solve this numerically... and very crude as well */
	cnt=0;
	rad = params->LJ_sig;
	step = -0.1 * params->LJ_sig;
	force=0.0;
	
	while(step != 0) {
	  frac2 = SQR(params->LJ_sig/rad);
	  frac6 = frac2*frac2*frac2;
	  force = 48.0 * params->LJ_eps * frac6*(frac6 - 0.5)/rad;
	  if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
	    step = - (step/2.0); 
	  }
	  if(fabs(force-force_cap) < 1.0e-6) step=0;
	  rad += step; cnt++;
	} 
      	params->LJ_capradius = rad;
      }
      else {
	params->LJ_capradius = 0.0; 
      }
      FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
			  this_node,i,j,rad,force,cnt));
    }
  }
}

#endif /* ifdef LENNARD_JONES */
