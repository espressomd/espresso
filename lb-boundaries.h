/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
 */

/** \file lb-boundaries.h
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.c.
 *
 */

#ifndef LB_BOUNDARIES_H
#define LB_BOUNDARIES_H

#include "utils.h"

#ifdef LB
#ifdef CONSTRAINTS

/** Initializes the constrains in the system. 
 *  This function determines the lattice sited which belong to boundaries
 *  and marks them with a corresponding flag. 
 */
void lb_init_constraints();

MDINLINE void lb_copy_neg_populations(LB_FluidNode *lbfluid, int to_index ,int from_index,double factor) {
	double weigth;
	int i;
#ifdef D3Q18
	weigth = -1.;
        for(i=0;i<18;i++)
          lbfluid[to_index].n[i]=factor*weigth*lbfluid[from_index].n[i];
#else
#error Boundary conditions are only implemented for D3Q18! (#defined in lb.h)
#endif
}

/** Apply boundary conditions to the LB fluid.
 * So far, only bounce-back boundary conditions are implemented.
 */
MDINLINE void lb_boundary_conditions() {

#ifdef D3Q18
  int k;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   1;                       // ( 1, 0, 0)
  next[1]  = - 1;                       // (-1, 0, 0)
  next[2]  =   yperiod;                 // ( 0, 1, 0)
  next[3]  = - yperiod;                 // ( 0,-1, 0)
  next[4]  =   zperiod;                 // ( 0, 0, 1)
  next[5]  = - zperiod;                 // ( 0, 0,-1)
  next[6]  =   (1+yperiod);             // ( 1, 1, 0)
  next[7]  = - (1+yperiod);             // (-1,-1, 0)
  next[8]  =   (1-yperiod);             // ( 1,-1, 0) 
  next[9]  = - (1-yperiod);             // (-1, 1, 0)
  next[10] =   (1+zperiod);             // ( 1, 0, 1)
  next[11] = - (1+zperiod);             // (-1, 0,-1)
  next[12] =   (1-zperiod);             // ( 1, 0,-1)
  next[13] = - (1-zperiod);             // (-1, 0, 1)
  next[14] =   (yperiod+zperiod);       // ( 0, 1, 1)
  next[15] = - (yperiod+zperiod);       // ( 0,-1,-1)
  next[16] =   (yperiod-zperiod);       // ( 0, 1,-1)
  next[17] = - (yperiod-zperiod);       // ( 0,-1, 1) 
  int reverse[] = { 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16 };

  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      /* bounce back to lower indices */
      lbfluid[k-next[0]].n[reverse[0]]   = lbfluid[k].n[0];
      lbfluid[k-next[2]].n[reverse[2]]   = lbfluid[k].n[2];
      lbfluid[k-next[4]].n[reverse[4]]   = lbfluid[k].n[4];
      lbfluid[k-next[6]].n[reverse[6]]   = lbfluid[k].n[6];
      lbfluid[k-next[9]].n[reverse[9]]   = lbfluid[k].n[9];
      lbfluid[k-next[10]].n[reverse[10]] = lbfluid[k].n[10];
      lbfluid[k-next[13]].n[reverse[13]] = lbfluid[k].n[13];
      lbfluid[k-next[14]].n[reverse[14]] = lbfluid[k].n[14];
      lbfluid[k-next[17]].n[reverse[17]] = lbfluid[k].n[17];

      lbfluid[k].n[0]  = 0.0;
      lbfluid[k].n[2]  = 0.0;
      lbfluid[k].n[4]  = 0.0;
      lbfluid[k].n[6]  = 0.0;
      lbfluid[k].n[9]  = 0.0;
      lbfluid[k].n[10] = 0.0;
      lbfluid[k].n[13] = 0.0;
      lbfluid[k].n[14] = 0.0;
      lbfluid[k].n[17] = 0.0;

    }

  }

  for (k=lblattice.halo_grid_volume-lblattice.halo_offset;k>=0;k--) {

    if (lbfluid[k].boundary) {

      /* bounce back to higher indices */
      lbfluid[k-next[1]].n[reverse[1]]   = lbfluid[k].n[1];
      lbfluid[k-next[3]].n[reverse[3]]   = lbfluid[k].n[3];
      lbfluid[k-next[5]].n[reverse[5]]   = lbfluid[k].n[5];
      lbfluid[k-next[7]].n[reverse[7]]   = lbfluid[k].n[7];
      lbfluid[k-next[8]].n[reverse[8]]   = lbfluid[k].n[8];
      lbfluid[k-next[11]].n[reverse[11]] = lbfluid[k].n[11];
      lbfluid[k-next[12]].n[reverse[12]] = lbfluid[k].n[12];
      lbfluid[k-next[15]].n[reverse[15]] = lbfluid[k].n[15];
      lbfluid[k-next[16]].n[reverse[16]] = lbfluid[k].n[16];

      lbfluid[k].n[1]  = 0.0;
      lbfluid[k].n[3]  = 0.0;
      lbfluid[k].n[5]  = 0.0;
      lbfluid[k].n[7]  = 0.0;
      lbfluid[k].n[8]  = 0.0;
      lbfluid[k].n[11] = 0.0;
      lbfluid[k].n[12] = 0.0;
      lbfluid[k].n[15] = 0.0;
      lbfluid[k].n[16] = 0.0;

    }

  }
#else
#error Boundary conditions are only implemented for D3Q18! (#defined in lb.h)
#endif
}

#endif /* CONSTRAINTS */
#endif /* LB */

#endif /* LB_BOUNDARIES_H */
