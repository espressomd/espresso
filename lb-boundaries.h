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
 * Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
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
#include "halo.h"

#ifdef LB
#ifdef CONSTRAINTS

#define LB_BOUNDARY_NONE                0
#define LB_BOUNDARY_BOUNCE_BACK         1
#define LB_BOUNDARY_SPECULAR_REFLECTION 2
#define LB_BOUNDARY_SLIP_REFLECTION     3
#define LB_BOUNDARY_PARTIAL_SLIP        4

typedef struct {
  int type;
  double slip_pref;
} LB_Boundary;

extern LB_Boundary lb_boundary_par;

/** Initializes the constrains in the system. 
 *  This function determines the lattice sited which belong to boundaries
 *  and marks them with a corresponding flag. 
 */
void lb_init_constraints();

/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 */
MDINLINE void lb_bounce_back() {

#ifdef D3Q19
  int k;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[19];
  next[0]  =   0;                       // ( 0, 0, 0) =
  next[1]  =   1;                       // ( 1, 0, 0) +
  next[2]  = - 1;                       // (-1, 0, 0)
  next[3]  =   yperiod;                 // ( 0, 1, 0) +
  next[4]  = - yperiod;                 // ( 0,-1, 0)
  next[5]  =   zperiod;                 // ( 0, 0, 1) +
  next[6]  = - zperiod;                 // ( 0, 0,-1)
  next[7]  =   (1+yperiod);             // ( 1, 1, 0) +
  next[8]  = - (1+yperiod);             // (-1,-1, 0)
  next[9]  =   (1-yperiod);             // ( 1,-1, 0) 
  next[10] = - (1-yperiod);             // (-1, 1, 0) +
  next[11] =   (1+zperiod);             // ( 1, 0, 1) +
  next[12] = - (1+zperiod);             // (-1, 0,-1)
  next[13] =   (1-zperiod);             // ( 1, 0,-1)
  next[14] = - (1-zperiod);             // (-1, 0, 1) +
  next[15] =   (yperiod+zperiod);       // ( 0, 1, 1) +
  next[16] = - (yperiod+zperiod);       // ( 0,-1,-1)
  next[17] =   (yperiod-zperiod);       // ( 0, 1,-1)
  next[18] = - (yperiod-zperiod);       // ( 0,-1, 1) +
  int reverse[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };

  /* bottom-up sweep */
  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      /* bounce back from lower indices */
      lbfluid[k].n[reverse[5]]  = lbfluid[k-next[5]].n[5];
      lbfluid[k].n[reverse[11]] = lbfluid[k-next[11]].n[11];
      lbfluid[k].n[reverse[14]] = lbfluid[k-next[14]].n[14];
      lbfluid[k].n[reverse[15]] = lbfluid[k-next[15]].n[15];
      lbfluid[k].n[reverse[18]] = lbfluid[k-next[18]].n[18];

      lbfluid[k-next[5]].n[5]   = 0.0;
      lbfluid[k-next[11]].n[11] = 0.0;
      lbfluid[k-next[14]].n[14] = 0.0;
      lbfluid[k-next[15]].n[15] = 0.0;
      lbfluid[k-next[18]].n[18] = 0.0;

    }

  }
  
  /* top-down sweep */
  for (k=lblattice.halo_grid_volume-lblattice.halo_offset;k>=0;k--) {

    if (lbfluid[k].boundary) {

      /* bounce back from higher indices */
      lbfluid[k].n[reverse[6]]  = lbfluid[k-next[6]].n[6];
      lbfluid[k].n[reverse[12]] = lbfluid[k-next[12]].n[12];
      lbfluid[k].n[reverse[13]] = lbfluid[k-next[13]].n[13];
      lbfluid[k].n[reverse[16]] = lbfluid[k-next[16]].n[16];
      lbfluid[k].n[reverse[17]] = lbfluid[k-next[17]].n[17];     

      lbfluid[k-next[6]].n[6]   = 0.0;   
      lbfluid[k-next[12]].n[12] = 0.0;  
      lbfluid[k-next[13]].n[13] = 0.0;  
      lbfluid[k-next[16]].n[16] = 0.0;  
      lbfluid[k-next[17]].n[17] = 0.0;

    }

  }

#else
#error Bounce back boundary conditions are only implemented for D3Q19!
#endif

}

/** Specular reflections at boundaries.
 * The populations that have propagated into a boundary node
 * are reflected in the direction perpendicular to the wall.
 * The tangential components of the populations are unchanged.
 * This results in full slip boundary conditions.
 */
MDINLINE void lb_specular_reflections() {

#ifdef D3Q19
  int k;
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[19];
  next[0]  =   0;                       // ( 0, 0, 0)
  next[1]  =   0;                       // ( 1, 0, 0)
  next[2]  = - 0;                       // (-1, 0, 0)
  next[3]  =   0;                       // ( 0, 1, 0)
  next[4]  = - 0;                       // ( 0,-1, 0)
  next[5]  =   zperiod;                 // ( 0, 0, 1)
  next[6]  = - zperiod;                 // ( 0, 0,-1)
  next[7]  =   0;                       // ( 1, 1, 0)
  next[8]  = - 0;                       // (-1,-1, 0)
  next[9]  =   0;                       // ( 1,-1, 0)
  next[10]  = - 0;                      // (-1, 1, 0)
  next[11] =   zperiod;                 // ( 1, 0, 1)
  next[12] = - zperiod;                 // (-1, 0,-1)
  next[13] = - zperiod;                 // ( 1, 0,-1)
  next[14] =   zperiod;                 // (-1, 0, 1)
  next[15] =   zperiod;                 // ( 0, 1, 1)
  next[16] = - zperiod;                 // ( 0,-1,-1)
  next[17] = - zperiod;                 // ( 0, 1,-1)
  next[18] =   zperiod;                 // ( 0,-1, 1)
  int reflect[] = { 0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 13, 14, 11, 12, 17, 18, 15, 16 };

  /* bottom-up sweep */
  for (k=zperiod;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      /* reflect from lower indices */
      lbfluid[k].n[reflect[5]]  = lbfluid[k-next[5]].n[5];
      lbfluid[k].n[reflect[11]] = lbfluid[k-next[11]].n[11];
      lbfluid[k].n[reflect[14]] = lbfluid[k-next[14]].n[14];
      lbfluid[k].n[reflect[15]] = lbfluid[k-next[15]].n[15];
      lbfluid[k].n[reflect[18]] = lbfluid[k-next[18]].n[18];

      lbfluid[k-next[5]].n[5]   = 0.0;
      lbfluid[k-next[11]].n[11] = 0.0;
      lbfluid[k-next[14]].n[14] = 0.0;
      lbfluid[k-next[15]].n[15] = 0.0;
      lbfluid[k-next[18]].n[18] = 0.0;

    }
    
  }

  /* top-down sweep */
  for (k=lblattice.halo_grid_volume-1-zperiod;k>=0;k--) {

    if (lbfluid[k].boundary) {

      /* reflect from higher indices */
      lbfluid[k].n[reflect[6]]  = lbfluid[k-next[6]].n[6];
      lbfluid[k].n[reflect[12]] = lbfluid[k-next[12]].n[12];
      lbfluid[k].n[reflect[13]] = lbfluid[k-next[13]].n[13];
      lbfluid[k].n[reflect[16]] = lbfluid[k-next[16]].n[16];
      lbfluid[k].n[reflect[17]] = lbfluid[k-next[17]].n[17];

      lbfluid[k-next[6]].n[6]   = 0.0;
      lbfluid[k-next[12]].n[12] = 0.0;
      lbfluid[k-next[13]].n[13] = 0.0;
      lbfluid[k-next[16]].n[16] = 0.0;
      lbfluid[k-next[17]].n[17] = 0.0;

    }
    
  }

#else
#error Specular reflections are only implemented for D3Q19!
#endif

}

/* slip reflection boundary condition */
MDINLINE void lb_slip_reflection() {
  int k;
  double s = lb_boundary_par.slip_pref;
  double r = 1.0 - s;

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   1;                           // ( 1, 0, 0) +
  next[1]  = - 1;                           // (-1, 0, 0)
  next[2]  =   yperiod;                     // ( 0, 1, 0) +
  next[3]  = - yperiod;                     // ( 0,-1, 0)
  next[4]  =   zperiod;                     // ( 0, 0, 1) +
  next[5]  = - zperiod;                     // ( 0, 0,-1)
  next[6]  =   (1+yperiod);                 // ( 1, 1, 0) +
  next[7]  = - (1+yperiod);                 // (-1,-1, 0)
  next[8]  =   (1-yperiod);                 // ( 1,-1, 0) 
  next[9]  = - (1-yperiod);                 // (-1, 1, 0) +
  next[10] =   (1+zperiod);                 // ( 1, 0, 1) +
  next[11] = - (1+zperiod);                 // (-1, 0,-1)
  next[12] =   (1-zperiod);                 // ( 1, 0,-1)
  next[13] = - (1-zperiod);                 // (-1, 0, 1) +
  next[14] =   (yperiod+zperiod);           // ( 0, 1, 1) +
  next[15] = - (yperiod+zperiod);           // ( 0,-1,-1)
  next[16] =   (yperiod-zperiod);           // ( 0, 1,-1)
  next[17] = - (yperiod-zperiod);           // ( 0,-1, 1) +

  int reflect[] = { 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 16, 17, 14, 15 };

  /* bottom up sweep */
  for (k=yperiod+zperiod; k<lblattice.halo_grid_volume; k++) {

    if (lbfluid[k].boundary == -1) {      

      lbfluid[k].n[reflect[4]]  = lbfluid[k-next[4]].n[4];
      lbfluid[k].n[reflect[10]] = s*lbfluid[k-zperiod].n[10] + r*lbfluid[k-next[13]].n[13];
      lbfluid[k].n[reflect[13]] = s*lbfluid[k-zperiod].n[13] + r*lbfluid[k-next[10]].n[10];
      lbfluid[k].n[reflect[14]] = s*lbfluid[k-zperiod].n[14] + r*lbfluid[k-next[17]].n[17];;
      lbfluid[k].n[reflect[17]] = s*lbfluid[k-zperiod].n[17] + r*lbfluid[k-next[14]].n[14];

      /* delete upgoing velocities (take care to delete the right ones!) */
      lbfluid[k-zperiod].n[4]  = 0.0;
      lbfluid[k-next[10]].n[10] = 0.0;
      lbfluid[k-zperiod].n[13] = 0.0;
      lbfluid[k-next[14]].n[14] = 0.0;
      lbfluid[k-zperiod].n[17] = 0.0;

      
    }

  }

  /* top down sweep */
  for (k=lblattice.halo_grid_volume-1-yperiod-zperiod; k>=0; k--) {

    if (lbfluid[k].boundary == 1) {

      lbfluid[k].n[reflect[5]]  = lbfluid[k-next[5]].n[5];
      lbfluid[k].n[reflect[11]] = s*lbfluid[k+zperiod].n[11] + r*lbfluid[k-next[12]].n[12];
      lbfluid[k].n[reflect[12]] = s*lbfluid[k+zperiod].n[12] + r*lbfluid[k-next[11]].n[11];
      lbfluid[k].n[reflect[15]] = s*lbfluid[k+zperiod].n[15] + r*lbfluid[k-next[16]].n[16];
      lbfluid[k].n[reflect[16]] = s*lbfluid[k+zperiod].n[16] + r*lbfluid[k-next[15]].n[15];

      /* delete downgoing velocities (take care to delete the right ones!) */
      lbfluid[k+zperiod].n[5]   = 0.0;
      lbfluid[k-next[11]].n[11] = 0.0;
      lbfluid[k+zperiod].n[12] = 0.0;
      lbfluid[k-next[15]].n[15] = 0.0;
      lbfluid[k+zperiod].n[16] = 0.0;
      
    }

  }

}

MDINLINE void lb_set_boundary_fields(LB_FluidNode *node, const double rho, const double *j, const double *pi) {

  switch (lb_boundary_par.type) {

  case LB_BOUNDARY_NONE:
  case LB_BOUNDARY_BOUNCE_BACK:
  case LB_BOUNDARY_SPECULAR_REFLECTION:
  case LB_BOUNDARY_SLIP_REFLECTION:
    lb_set_local_fields(node, 0.0, j, pi);
    break;
  case LB_BOUNDARY_PARTIAL_SLIP:
    lb_boundary_equilibrium2(node, rho, j, pi);
    break;

  }

}

/** Apply boundary conditions to the LB fluid. */
MDINLINE void lb_boundary_conditions() {

  switch (lb_boundary_par.type) {

  case LB_BOUNDARY_NONE:
    break;
  case LB_BOUNDARY_BOUNCE_BACK:
    lb_bounce_back();
    break;
  case LB_BOUNDARY_SPECULAR_REFLECTION:
    lb_specular_reflections();
    break;
  case LB_BOUNDARY_SLIP_REFLECTION:
    lb_slip_reflection();
    break;
  case LB_BOUNDARY_PARTIAL_SLIP:
    lb_boundary_equilibration();
    break;

  }

}

/** Parser for the \ref lbfluid command. */
int lbboundaries_cmd(ClientData data, Tcl_Interp *interp, int argc, char **argv);

#endif /* CONSTRAINTS */
#endif /* LB */

#endif /* LB_BOUNDARIES_H */
