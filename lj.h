// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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

MDINLINE void add_lj_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  int j;
  double r_off, frac2, frac6, fac=0.0;
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) { 
    r_off = dist - ia_params->LJ_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);

      for(j=0;j<3;j++) {
	p1->f.f[j] += fac * d[j];
	p2->f.f[j] -= fac * d[j];
      }
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
      for(j=0;j<3;j++) {
	/* vector d is rescaled to length LJ_capradius */
	p1->f.f[j] += fac * d[j];
	p2->f.f[j] -= fac * d[j];
      }
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / ia_params->LJ_capradius;

      p1->f.f[0] += fac * ia_params->LJ_capradius;
      p2->f.f[0] -= fac * ia_params->LJ_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}

MDINLINE double lj_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;

  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
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

#endif
