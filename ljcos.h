// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef LJCOS_H
#define LJCOS_H

/** \file ljcos.h
 *  Routines to calculate the lennard jones+cosine energy and/or force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
*/

MDINLINE void add_ljcos_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  int j;
  double r_off, frac2, frac6, fac=0.0;

  if(dist < ia_params->LJCOS_cut+ia_params->LJCOS_offset) {
    r_off = dist - ia_params->LJCOS_offset;
    /* cos part of ljcos potential. */
    if(dist > ia_params->LJCOS_rmin+ia_params->LJCOS_offset) {
      fac   = (r_off/dist) * ia_params->LJCOS_alfa * ia_params->LJCOS_eps * (sin(ia_params->LJCOS_alfa * SQR(r_off) + ia_params->LJCOS_beta));
      for(j=0;j<3;j++) {
	    /* vector d is rescaled to length LJ_capradius */
	    p1->f[j] += fac * d[j];
	    p2->f[j] -= fac * d[j];
      }
    }
    /* lennard-jones part of the potential. */
    else if(dist > 0) {
      frac2 = SQR(ia_params->LJCOS_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJCOS_eps * frac6*(frac6 - 0.5) / (r_off * dist);

      for(j=0;j<3;j++) {
	    p1->f[j] += fac * d[j];
	    p2->f[j] -= fac * d[j];
      }
#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJCOS-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->r.identity,p2->r.identity,fac*dist,dist);
#endif
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->r.identity,p2->r.identity));

      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / ia_params->LJ_capradius;

      p1->f[0] += fac * ia_params->LJ_capradius;
      p2->f[0] -= fac * ia_params->LJ_capradius;
    }

    ONEPART_TRACE(if(p1->r.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f[0],p1->f[1],p1->f[2],p2->r.identity,dist,fac));
    ONEPART_TRACE(if(p2->r.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f[0],p2->f[1],p2->f[2],p1->r.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->r.identity,p2->r.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}


MDINLINE double ljcos_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;

  if(dist < ia_params->LJCOS_cut+ia_params->LJCOS_offset) {
    r_off = dist-ia_params->LJCOS_offset;
    /* lennard-jones part of the potential. */
    if (dist < (ia_params->LJCOS_rmin+ia_params->LJCOS_offset)) {
      //printf("this is nomal ,  %.3e \n",r_off);
      frac2 = SQR(ia_params->LJCOS_sig/r_off);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJCOS_eps*(SQR(frac6)-frac6);
    }
    /* cosine part of the potential. */
    else if (dist < (ia_params->LJCOS_cut+ia_params->LJCOS_offset)) {
      return .5*ia_params->LJCOS_eps*(cos(ia_params->LJCOS_alfa*SQR(r_off)+ia_params->LJCOS_beta)-1.);
    }
    /* this should not happen! */
    else {
      fprintf(stderr,"this is the distance, which is negative %.3e\n",r_off);
    }
  }
  return 0.0;
}


#endif
