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

MDINLINE void add_fene_pair_force(Particle *p1, Particle *p2, int type_num)
{
  int i;
  double dx[3], dist2=0.0, fac;
  for(i=0;i<3;i++) {
    dx[i] = p1->r.p[i] - p2->r.p[i];
    dx[i] -= dround(dx[i]/box_l[i])*box_l[i];
    dist2 += SQR(dx[i]);
  }
  
  FENE_TRACE(if(dist2 >= SQR(bonded_ia_params[type_num].p.fene.r_fene)) fprintf(stderr,"FENE Bond between Pair (%d,%d) broken: dist=%f\n",p1->identity,p2->identity,sqrt(dist2)) );
  
  fac = bonded_ia_params[type_num].p.fene.k_fene;
  fac /= (1.0 - dist2/SQR(bonded_ia_params[type_num].p.fene.r_fene));

  FENE_TRACE(if(fac > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", p1->identity,p2->identity,fac,sqrt(dist2)) );

  for(i=0;i<3;i++) {
    p1->f[i] -= fac*dx[i];
    p2->f[i] += fac*dx[i];
  }
}

#endif
