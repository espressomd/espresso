/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
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
