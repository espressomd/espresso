// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef SUBT_LJ_H
#define SUBT_LJ_H
/** \file subt_lj.h
 *  Routines to subtract the LENNARD-JONES Energy and/or the LENNARD-JONES force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sathish@mpip-mainz.mpg.de">sathish</a>
*/

/************************************************************/

/** Computes the negative of the LENNARD-JONES pair forces 
    and adds this force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param type_num  bond type number of this interaction (see \ref #inter).
*/
MDINLINE void add_subt_lj_pair_force(Particle *p1, Particle *p2, int type_num)
{
  int i;
  double dx[3], dist=0.0, dist2=0.0, fac_lj=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;

  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  dist=sqrt(dist2);

  if(dist >= bonded_ia_params[type_num].p.subt_lj.r) {
    fprintf(stderr,"%d: add_lj_pair_force: ERROR: Bond between Pair (%d,%d) is too long: dist=%f\n",this_node,
	    p1->p.identity,p2->p.identity,dist); 
    errexit();
  }

  ia_params = get_ia_param(p1->p.type,p2->p.type);
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

  
  for(i=0;i<3;i++) {
    p1->f.f[i] -= fac_lj*dx[i];
    p2->f.f[i] += fac_lj*dx[i];
#ifdef NPT
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir -= fac_lj*dx[i] * dx[i];
#endif
  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac_lj %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(dist2),fac_lj));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac_lj %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(dist2),fac_lj));
}

MDINLINE double subt_lj_pair_energy(Particle *p1, Particle *p2, int type_num)
{
  double dx[3], dist2=0.0, dist=0.0, energy_lj=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  dist=sqrt(dist2);  
  
  if(dist >= bonded_ia_params[type_num].p.subt_lj.r) {
    fprintf(stderr,"%d: add_subt_lj_pair_energy: ERROR: Bond between Pair (%d,%d) is too long: dist=%f\n",
	    this_node,p1->p.identity,p2->p.identity,dist); 
    errexit();
    	}
  
  ia_params = get_ia_param(p1->p.type,p2->p.type);
  
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
    return -energy_lj;
  }

#endif
