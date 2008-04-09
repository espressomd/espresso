// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef LJ_H
#define LJ_H

/** \file lj.h
 *  Routines to calculate the lennard jones energy and/or  force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/

#ifdef LENNARD_JONES

MDINLINE int printljIAToResult(Tcl_Interp *interp, int i, int j)
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
  
  return TCL_OK;
}

/** set the force cap for the LJ interaction.
    @param ljforcecap the maximal force, 0 to disable, -1 for individual cutoff
    for each of the interactions.
*/
MDINLINE int ljforcecap_set_params(double ljforcecap)
{
  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);
  
  return TCL_OK;
}

MDINLINE int lennard_jones_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double cut,
				      double shift, double offset,
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

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}

/// parser for the forcecap
MDINLINE int inter_parse_ljforcecap(Tcl_Interp * interp, int argc, char ** argv)
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

MDINLINE int lj_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJ */
  double eps, sig, cut, shift, offset, cap_radius;
  int change;

  /* get lennard-jones interaction type */
  if (argc < 6) {
    Tcl_AppendResult(interp, "lennard-jones needs 5 parameters: "
		     "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
		     (char *) NULL);
    return 0;
  }

  /* copy lennard-jones parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, sig))   ||
      (! ARG_IS_D(3, cut))   ||
      (! ARG_IS_D(4, shift)) ||
      (! ARG_IS_D(5, offset)    )) {
    Tcl_AppendResult(interp, "lennard-jones needs 5 DOUBLE parameters: "
		     "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 6;
	
  cap_radius = -1.0;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 7 && ARG_IS_D(6, cap_radius))
    change++;
  else
    Tcl_ResetResult(interp);
  if (lennard_jones_set_params(part_type_a, part_type_b,
			       eps, sig, cut, shift, offset,
			       cap_radius) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}


/** Calculate lennard Jones force between particle p1 and p2 */
MDINLINE void add_lj_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double r_off, frac2, frac6, fac=0.0;

#ifdef WATER
  double p1_com[3],p2_com[3],com_dist;

  if (p1->p.mol_id==p2->p.mol_id) return;

  if ((get_com_h2o(p1,p1_com) == -1 ) || (get_com_h2o(p2,p2_com)==-1)){
     return;
  }
  else{
     com_dist=min_distance(p1_com,p2_com);
  }

  if (com_dist < ia_params->LJ_cut+ia_params->LJ_offset){
#else
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
#endif
    r_off = dist - ia_params->LJ_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);

      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJ-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
      for(j=0;j<3;j++)
	/* vector d is rescaled to length LJ_capradius */
	force[j] += fac * d[j];
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / ia_params->LJ_capradius;

      force[0] += fac * ia_params->LJ_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}

/** calculate Lennard jones energy between particle p1 and p2. */
MDINLINE double lj_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;

#ifdef WATER
  double p1_com[3],p2_com[3],com_dist;

  if (p1->p.mol_id==p2->p.mol_id) return 0.0;

  if ((get_com_h2o(p1,p1_com) == -1 ) || (get_com_h2o(p2,p2_com)==-1)){
     return 0.0;
  }
  else{
     com_dist=min_distance(p1_com,p2_com);
  }

  if (com_dist < ia_params->LJ_cut+ia_params->LJ_offset){
#else
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
#endif
    r_off = dist - ia_params->LJ_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
    /* this should not happen! */
    else {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
  }
  return 0.0;
}

/** calculate lj_capradius from lj_force_cap */
MDINLINE void calc_lj_cap_radii(double force_cap)
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
#endif
