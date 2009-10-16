// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2007; all rights reserved unless otherwise stated.

/** \file ljgen.c
 *  Routines to calculate the generalized lennard jones energy and/or  force 
 *  for a particle pair. "Generalized" here means that the LJ energy is of the
 *  form
 *
 *  eps * [ b1 * (sigma/(r-r_offset))^a1 - b2 * (sigma/(r-r_offset))^a2 + shift]
 *
 *  \ref forces.c
*/

#include "config.h"

#ifdef LENNARD_JONES_GENERIC

#include "ljgen.h"
#include "communication.h"
#include "parser.h"

int printljgenIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->LJGEN_eps, buffer);
  Tcl_AppendResult(interp, "lj-gen ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_shift, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  snprintf (buffer, sizeof (buffer), "%d ", data->LJGEN_a1);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  snprintf (buffer, sizeof (buffer), "%d ", data->LJGEN_a2);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  Tcl_PrintDouble(interp, data->LJGEN_b1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_b2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_capradius, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
 
  return TCL_OK;
}


static int ljgen_set_params(int part_type_a, int part_type_b,
			       double eps, double sig, double cut,
			       double shift, double offset,
			       int a1, int a2, double b1, double b2,
			       double cap_radius)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJGEN should be symmetrically */
  data->LJGEN_eps    = data_sym->LJGEN_eps    = eps;
  data->LJGEN_sig    = data_sym->LJGEN_sig    = sig;
  data->LJGEN_cut    = data_sym->LJGEN_cut    = cut;
  data->LJGEN_shift  = data_sym->LJGEN_shift  = shift;
  data->LJGEN_offset = data_sym->LJGEN_offset = offset;
  data->LJGEN_a1     = data_sym->LJGEN_a1 = a1;
  data->LJGEN_a2     = data_sym->LJGEN_a2 = a2;
  data->LJGEN_b1     = data_sym->LJGEN_b1 = b1;
  data->LJGEN_b2     = data_sym->LJGEN_b2 = b2;
 
  if (cap_radius > 0) {
    data->LJGEN_capradius = cap_radius;
    data_sym->LJGEN_capradius = cap_radius;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}


int ljgen_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJGEN */
  double eps, sig, cut, shift, offset, cap_radius, b1, b2;
  int change, a1, a2;

  /* get lennard-jones interaction type */
  if (argc < 10) {
    Tcl_AppendResult(interp, "lj-gen needs 9 parameters: "
		     "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset> <a1> <a2> <b1> <b2>",
		     (char *) NULL);
    return 0;
  }

  /* copy lennard-jones parameters */
  if ((! ARG_IS_D(1, eps))    ||
      (! ARG_IS_D(2, sig))    ||
      (! ARG_IS_D(3, cut))    ||
      (! ARG_IS_D(4, shift))  ||
      (! ARG_IS_D(5, offset)) ||
      (! ARG_IS_I(6, a1))     ||
      (! ARG_IS_I(7, a2))     ||
      (! ARG_IS_D(8, b1))     ||
      (! ARG_IS_D(9, b2))) {
    Tcl_AppendResult(interp, "lj-gen needs 7 DOUBLE and 2 INT parameers: "
		     "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset> <a1> <a2> <b1> <b2>",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 10;
	
  cap_radius = -1.0;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 11 && ARG_IS_D(10, cap_radius))
    change++;
  else
    Tcl_ResetResult(interp);
  if (ljgen_set_params(part_type_a, part_type_b,
		       eps, sig, cut, shift, offset, a1, a2, b1, b2,
		       cap_radius) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}


/** double**integer power function.  */
MDINLINE double powi (double x, int n)
{
  double y = n % 2 ? x : 1;
  while (n >>= 1)
    {
      x = x * x;
      if (n % 2)
	y = y * x;
    }
  return y;
}


/** Calculate lennard Jones force between particle p1 and p2 */
void add_ljgen_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double r_off, frac, fac=0.0;
  if(dist < ia_params->LJGEN_cut+ia_params->LJGEN_offset) { 
    r_off = dist - ia_params->LJGEN_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJGEN_capradius) {
      frac = ia_params->LJGEN_sig/r_off;
      fac   = ia_params->LJGEN_eps 
	* (ia_params->LJGEN_b1 * ia_params->LJGEN_a1 * powi(frac, ia_params->LJGEN_a1) 
	   - ia_params->LJGEN_b2 * ia_params->LJGEN_a2 * powi(frac, ia_params->LJGEN_a2))
	/ (r_off * dist);
      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJ-Gen-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      fac   = ia_params->LJGEN_eps 
	* (ia_params->LJGEN_b1 * ia_params->LJGEN_a1 * powi(frac, ia_params->LJGEN_a1) 
	   - ia_params->LJGEN_b2 * ia_params->LJGEN_a2 * powi(frac, ia_params->LJGEN_a2))
	/ (ia_params->LJGEN_capradius * dist);
      for(j=0;j<3;j++)
	/* vector d is rescaled to length LJGEN_capradius */
	force[j] += fac * d[j];
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones-Generic warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      fac   = ia_params->LJGEN_eps 
	* (ia_params->LJGEN_b1 * ia_params->LJGEN_a1 * powi(frac, ia_params->LJGEN_a1) 
	   - ia_params->LJGEN_b2 * ia_params->LJGEN_a2 * powi(frac, ia_params->LJGEN_a2))
	/ ia_params->LJGEN_capradius;

      force[0] += fac * ia_params->LJGEN_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJGEN   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJGEN   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJGEN: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}

/** calculate Lennard jones energy between particle p1 and p2. */
double ljgen_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac;

  if(dist < ia_params->LJGEN_cut+ia_params->LJGEN_offset) {
    r_off = dist - ia_params->LJGEN_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJGEN_capradius) {
      frac = ia_params->LJGEN_sig/r_off;
      return ia_params->LJGEN_eps*(
	      ia_params->LJGEN_b1 * powi(frac, ia_params->LJGEN_a1) 
	      - ia_params->LJGEN_b2 * powi(frac, ia_params->LJGEN_a2)
				    + ia_params->LJGEN_shift);
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      return ia_params->LJGEN_eps*(
	      ia_params->LJGEN_b1 * powi(frac, ia_params->LJGEN_a1) 
	      - ia_params->LJGEN_b2 * powi(frac, ia_params->LJGEN_a2)
				    + ia_params->LJGEN_shift);
    }
    /* this should not happen! */
    else {
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      return ia_params->LJGEN_eps*(
	      ia_params->LJGEN_b1 * powi(frac, ia_params->LJGEN_a1) 
	       - ia_params->LJGEN_b2 * powi(frac, ia_params->LJGEN_a2)
				    + ia_params->LJGEN_shift);
    }
  }
  return 0.0;
}

/** calculate lj_capradius from lj_force_cap */
void calc_ljgen_cap_radii(double force_cap)
{
  int i,j,cnt=0;
  IA_parameters *params;
  double force=0.0, rad=0.0, step, frac;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      if(force_cap > 0.0 && params->LJGEN_eps > 0.0) {
	/* I think we have to solve this numerically... and very crude as well */
	cnt=0;
	rad = params->LJGEN_sig;
	step = -0.1 * params->LJGEN_sig;
	force=0.0;
	
	while(step != 0) {
	  frac = params->LJGEN_sig/rad;
	  force =  params->LJGEN_eps 
	    * (params->LJGEN_b1 * params->LJGEN_a1 * powi(frac, params->LJGEN_a1) 
	       - params->LJGEN_b2 * params->LJGEN_a2 * powi(frac, params->LJGEN_a2))/rad;
	  if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
	    step = - (step/2.0); 
	  }
	  if(fabs(force-force_cap) < 1.0e-6) step=0;
	  rad += step; cnt++;
	} 
      	params->LJGEN_capradius = rad;
      }
      else {
	params->LJGEN_capradius = 0.0; 
      }
      FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
			  this_node,i,j,rad,force,cnt));
    }
  }
}

#endif
