// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef SUBT_LJ_HARM_H
#define SUBT_LJ_HARM_H
/** \file SUBT_LJ_HARM.h
 *  Routines to subtract the LENNARD-JONES Energy and/or the LENNARD-JONES force 
 *  from the HARMONIC Energy and/or the force for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sathish@mpip-mainz.mpg.de">sathish</a>
*/

/************************************************************/

/** Computes the difference between the HARMONIC and the LENNARD-JONES pair forces 
    and adds this force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param type_num  bond type number of this interaction (see \ref #inter).
*/
MDINLINE void add_subt_lj_harm_pair_force(Particle *p1, Particle *p2, int type_num)
{
  int i;
  double dx[3], dist=0.0, dist2=0.0, fac_harm=0.0, fac_lj=0.0, fac;
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  dist=sqrt(dist2);

  fac_harm = bonded_ia_params[type_num].p.subt_lj_harm.k;
  fac_harm *= (dist-bonded_ia_params[type_num].p.subt_lj_harm.r);
  fac_harm /= dist;

  double r_off, frac2, frac6;
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) { 
    r_off = dist - ia_params->LJ_offset;

    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac_lj   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);			  
    	}

    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac_lj   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
      	    }      
    }	    

    fac = fac_harm + fac_lj;
    	
    for(i=0;i<3;i++) {
    	    p1->f.f[i] -= fac*dx[i];
    	    p2->f.f[i] += fac*dx[i];
    	    }

  ONEPART_TRACE(if(p1->r.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ_HARM f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f[0],p1->f[1],p1->f[2],p2->r.identity,sqrt(dist2),fac));
  ONEPART_TRACE(if(p2->r.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ_HARM f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f[0],p2->f[1],p2->f[2],p1->r.identity,sqrt(dist2),fac));
}

MDINLINE double subt_lj_harm_pair_energy(Particle *p1, Particle *p2, int type_num)
{
  double dx[3], dist2=0.0, dist=0.0, energy_harm=0.0, energy_lj=0.0;
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  dist=sqrt(dist2);
  
  energy_harm = 0.5*bonded_ia_params[type_num].p.subt_lj_harm.k;
  energy_harm *= SQR(dist-bonded_ia_params[type_num].p.subt_lj_harm.r);
  
  double r_off, frac2, frac6;

  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
    r_off = dist - ia_params->LJ_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      energy_lj = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    	    }   
     /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      energy_lj = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
      	    }
    	}
    return energy_harm-energy_lj;
  }

#endif
