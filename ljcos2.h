// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#ifndef LJCOS2_H
#define LJCOS2_H

/** \file ljcos2.h
 *  Routines to calculate the lennard-jones with cosine tail energy and/or  force 
 *  for a particle pair.  Cosine tail is different from that in ljcos.h
 *  Used for attractive tail/tail interactions in lipid bilayer calculations
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:reynolds@mpip-mainz.mpg.de">Ben</a>
*/

#ifdef LJCOS2
#include <math.h>

MDINLINE int printljcos2IAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->LJCOS2_eps, buffer);
  Tcl_AppendResult(interp, "lj-cos2 ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJCOS2_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJCOS2_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJCOS2_w, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  
  return TCL_OK;
}

MDINLINE int ljcos2_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double offset,
				      double w)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* lj-cos2 should be symmetrically */
  data->LJCOS2_eps    = data_sym->LJCOS2_eps    = eps;
  data->LJCOS2_sig    = data_sym->LJCOS2_sig    = sig;
  data->LJCOS2_offset = data_sym->LJCOS2_offset = offset;
  data->LJCOS2_w      = data_sym->LJCOS2_w      = w;

  /* calculate dependent parameters */
  data->LJCOS2_rchange = data_sym->LJCOS2_rchange = pow(2,1/6.)*sig;
  data->LJCOS2_cut     = data_sym->LJCOS2_cut     = w + data_sym->LJCOS2_rchange;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}


MDINLINE int ljcos2_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for lj-cos2 */
  double eps, sig, offset, w, cap_radius;
  int change;

  /* get lj-cos2 interaction type */
  if (argc < 5) {
    Tcl_AppendResult(interp, "ljcos2 needs 4 parameters: "
		     "<ljcos2_eps> <ljcos2_sig> <ljcos2_offset> <ljcos2_w>",
		     (char *) NULL);
    return 0;
  }

  /* copy lj-cos2 parameters */
  if ((! ARG_IS_D(1, eps))    ||
      (! ARG_IS_D(2, sig))    ||
      (! ARG_IS_D(3, offset)) ||
      (! ARG_IS_D(4, w))) {
    Tcl_AppendResult(interp, "ljcos2 needs 4 DOUBLE parameters: "
		     "<ljcos2_eps> <ljcos2_sig> <ljcos2_offset> <ljcos2_w>",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 5;

  cap_radius = -1;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 6 && ARG_IS_D(5, cap_radius))
    change++;
  else
    Tcl_ResetResult(interp);

  if (ljcos2_set_params(part_type_a, part_type_b,
			       eps, sig, offset, w
			       ) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}


/** Calculate lj-cos2 force between particle p1 and p2 */
MDINLINE void add_ljcos2_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double r_off, frac2, frac6, fac=0.0;
  if(dist < ia_params->LJCOS2_cut+ia_params->LJCOS2_offset) { 
    r_off = dist - ia_params->LJCOS2_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJCOS2_capradius) {
      if(r_off < ia_params->LJCOS2_rchange) {
        frac2 = SQR(ia_params->LJCOS2_sig/r_off);
        frac6 = frac2*frac2*frac2;
        fac   = 48.0 * ia_params->LJCOS2_eps * frac6*(frac6 - 0.5) / (r_off*dist);
      }
      else if (r_off< ia_params->LJCOS2_rchange + ia_params->LJCOS2_w) {
        fac   = -ia_params->LJCOS2_eps*M_PI/2/ia_params->LJCOS2_w/dist * sin(M_PI*(r_off-ia_params->LJCOS2_rchange)/ia_params->LJCOS2_w);
      }
      
      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJ-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    /* capped part of lj-cos2 potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJCOS2_sig/ia_params->LJCOS2_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJCOS2_eps * frac6*(frac6 - 0.5) / (ia_params->LJCOS2_capradius * dist);
      for(j=0;j<3;j++)
	/* vector d is rescaled to length LJCOS2_capradius */
	force[j] += fac * d[j];
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      frac2 = SQR(ia_params->LJCOS2_sig/ia_params->LJCOS2_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJCOS2_eps * frac6*(frac6 - 0.5) / ia_params->LJCOS2_capradius;

      force[0] += fac * ia_params->LJCOS2_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}

/** calculate lj-cos2 energy between particle p1 and p2. */
MDINLINE double ljcos2_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;

  if(dist < ia_params->LJCOS2_cut+ia_params->LJCOS2_offset) {
    r_off = dist - ia_params->LJCOS2_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJCOS2_capradius) {
      if (r_off < ia_params->LJCOS2_offset + ia_params->LJCOS2_rchange){
        frac2 = SQR(ia_params->LJCOS2_sig/r_off);
        frac6 = frac2*frac2*frac2;
        return 4.0*ia_params->LJCOS2_eps*(SQR(frac6)-frac6);
      }
      else if (r_off < ia_params->LJCOS2_rchange + ia_params->LJCOS2_w){
        return -ia_params->LJCOS2_eps/2 * (cos(M_PI*(r_off-ia_params->LJCOS2_rchange)/ia_params->LJCOS2_w)+1);
      }
    }
    /* capped part of lj-cos2 potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJCOS2_sig/ia_params->LJCOS2_capradius);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJCOS2_eps*(SQR(frac6)-frac6);
    }
    /* this should not happen! */
    else {
      frac2 = SQR(ia_params->LJCOS2_sig/ia_params->LJCOS2_capradius);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJCOS2_eps*(SQR(frac6)-frac6);
    }
  }
  return 0.0;
}

/** calculate ljcos2_capradius from ljcos2_force_cap */
MDINLINE void calc_ljcos2_cap_radii(double force_cap)
{
  int i,j,cnt=0;
  IA_parameters *params;
  double force=0.0, rad=0.0, step, frac2, frac6;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      if(force_cap > 0.0 && params->LJCOS2_eps > 0.0) {
	/* I think we have to solve this numerically... and very crude as well */
	cnt=0;
	rad = params->LJCOS2_sig;
	step = -0.1 * params->LJCOS2_sig;
	force=0.0;
	
	while(step != 0) {
	  frac2 = SQR(params->LJCOS2_sig/rad);
	  frac6 = frac2*frac2*frac2;
          if (rad < params->LJCOS2_rchange) {
            force = 48.0 * params->LJCOS2_eps * frac6*(frac6 - 0.5)/rad;
          }
          else {
	    force = 0;
	  }
	  if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
	    step = - (step/2.0); 
	  }
	  if(fabs(force-force_cap) < 1.0e-6) step=0;
	  rad += step; cnt++;
	} 
      	params->LJCOS2_capradius = rad;
      }
      else {
	params->LJCOS2_capradius = 0.0; 
      }
      FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
			  this_node,i,j,rad,force,cnt));
    }
  }
}

#endif /* ifdef LJCOS2 */
#endif
