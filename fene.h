#ifndef FENE_H
#define FENE_H
/** \file fene.h
 *  Routines to calculate the FENE Energy or/and FENE force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "communication.h"
#include "interaction_data.h"
#include "grid.h"

/************************************************************/

/* DEBUG FLAG */
#define FENE_DEBUG

#ifdef FENE_DEBUG
#define FENE_TRACE(cmd) { cmd; }
#else
#define FENE_TRACE(cmd)
#endif

/************************************************************/

MDINLINE void add_fene_pair_force(int p1_ind, int p2_ind, int type_num)
{
  int i;
  double dx[3], dist2=0.0, fac;
  if(p2_ind == -1) {
    fprintf(stderr,"%d: ERROR: Bonded atoms %d and %d not on the same node\n"
	    ,this_node, particles[p1_ind].identity,particles[p1_ind].bonds[i+1]); 
    errexit();
  }
  for(i=0;i<3;i++) {
    dx[i] = particles[p1_ind].p[i] - particles[p2_ind].p[i];
    dx[i] -= dround(dx[i]/box_l[i])*box_l[i];
    dist2 += SQR(dx[i]);
  }

  FENE_TRACE(if(dist2 >= SQR(bonded_ia_params[type_num].p.fene.r_fene)) fprintf(stderr,"FENE Bond between Pair (%d,%d) broken: dist=%f\n",particles[p1_ind].identity,particles[p2_ind].identity,sqrt(dist2)) );
  
  fac = bonded_ia_params[type_num].p.fene.k_fene;
  fac /= (1.0 - dist2/SQR(bonded_ia_params[type_num].p.fene.r_fene));

  FENE_TRACE(if(fac > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", particles[p1_ind].identity,particles[p2_ind].identity,fac,sqrt(dist2)) );

  for(i=0;i<3;i++) {
    particles[p1_ind].f[i] -= fac*dx[i];
    particles[p2_ind].f[i] += fac*dx[i];
  }
}

#endif
