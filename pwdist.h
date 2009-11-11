// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted 
// upon receiving the distribution and by which you are legally bound while 
// utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or 
// FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current 
// version can be found, or write to Max-Planck-Institute for Polymer Research, 
// Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

/*
#ifndef ANGLEDIST_H
#define ANGLEDIST_H
#ifndef ENDANGLEDIST_H
#define ENDANGLEDIST_H
*/

#ifdef BOND_ANGLEDIST
#ifdef BOND_ENDANGLEDIST

#include "utils.h"

/************************************************************/

/* Calculates the minimum distance between a particle to any wall constraint
 */
MDINLINE double calc_pwdist(Particle *p1, Bonded_ia_parameters *iaparams)
{
  int j,k,img[3];
  double distwallmin, distmx;
  double folded_pos_p1[3], vec[3];
  double pwdist[n_constraints];

#ifdef BOND_ENDANGLEDIST
  distmx = iaparams->p.endangledist.distmax;
#endif
#ifdef BOND_ANGLEDIST
  distmx = iaparams->p.angledist.distmax;
#endif

  /* folds coordinates of p_left into original box */
  memcpy(folded_pos_p2, p1->r.p, 3*sizeof(double));
  memcpy(img, p1->l.i, 3*sizeof(int));
  fold_position(folded_pos_p1, img);

  /* Gets and tests wall data */
  for(k=0;k<n_constraints;k++) {
    switch(constraints[k].type) {
      case CONSTRAINT_WAL: 
      /* wall normal */
      wall=constraints[k].c.wal;
      /* fprintf(stdout,"Wall normal %d = (%f,%f,%f)\n",k,wall.n[0],wall.n[1],wall.n[2]); */
      /* check that constraint vector is normalised */
      for(j=0;j<3;j++) normal += wall.n[j] * wall.n[j];
      if (sqrt(normal) != 1.0) {
        for(j=0;j<3;j++) wall.n[j]=wall.n[j]/normal;
      }
      break;
    }
  }

  /* Calculate distance of end particle from closest wall */
  for(k=0;k<n_constraints;k++) {
    switch(constraints[k].type) {
      case CONSTRAINT_WAL: 
      /* distwallmin is distance of closest wall from p1 */
      pwdist[k]=-1.0 * constraints[k].c.wal.d;
      for(j=0;j<3;j++) {
        pwdist[k] += folded_pos_end[j] * constraints[k].c.wal.n[j];
      }
      if (k==0) {
        distwallmin=pwdist[0];
      } else {
        if (pwdist[k] < distwallmin) {
          distwallmin = pwdist[k];
        }
      }
      /* fprintf(stdout,"Part[%d]-wall[%d] distance=%f\n",p1->p.identity,k,pwdist[k]); */
      break;
    }
  }

  return distwallmin;
}

#endif
#endif
