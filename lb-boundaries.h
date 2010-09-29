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
#include <tcl.h>
#include "utils.h"
#include "halo.h"
#include "constraint.h"
#include "config.h"
#include "lb.h"

#ifdef LB_BOUNDARIES

/** No constraint applied */
#define LB_BOUNDARY_NONE                0
#define LB_BOUNDARY_BOUNCE_BACK         1
#define LB_BOUNDARY_SPECULAR_REFLECTION 2
#define LB_BOUNDARY_SLIP_REFLECTION     3
#define LB_BOUNDARY_PARTIAL_SLIP        4

/** wall constraint applied */
#define LB_BOUNDARY_WAL 1
/** spherical constraint applied */
#define LB_BOUNDARY_SPH 2
/** (finite) cylinder shaped constraint applied */
#define LB_BOUNDARY_CYL 3

/** Structure to specify a boundary. */
typedef struct {
  /** type of the boundary. */
  int type;
  double slip_pref;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
  } c;
} LB_Boundary;
/*@}*/

//extern LB_Boundary lb_boundary_par;

MDINLINE void lb_calc_modes();

/** Initializes the constrains in the system. 
 *  This function determines the lattice sited which belong to boundaries
 *  and marks them with a corresponding flag. 
 */
void lb_init_boundaries();
int lb_boundary(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv);

/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 */
MDINLINE void lb_bounce_back() {

#ifdef D3Q19
#ifndef PULL
<<<<<<< lb-boundaries.h
  int k;
=======
  int k,i;
>>>>>>> 2.3.6.7
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

    if (lbfields[k].boundary) {
      for (i=0; i<19; i++) {
        lbfluid[1][reverse[i]][k-next[i]]   = lbfluid[1][i][k];
      }

//      /* bounce back to lower indices */
//      lbfluid[1][reverse[5]][k-next[5]]   = lbfluid[1][5][k];
//      lbfluid[1][reverse[11]][k-next[11]] = lbfluid[1][11][k];
//      lbfluid[1][reverse[14]][k-next[14]] = lbfluid[1][14][k];
//      lbfluid[1][reverse[15]][k-next[15]] = lbfluid[1][15][k];
//      lbfluid[1][reverse[18]][k-next[18]] = lbfluid[1][18][k];
//
//      /* delete populations in the wall */
//      lbfluid[1][5][k]  = - lbmodel.coeff[5][0]*lbpar.rho_lb_units; 
//      lbfluid[1][11][k] = - lbmodel.coeff[11][0]*lbpar.rho_lb_units;
//      lbfluid[1][14][k] = - lbmodel.coeff[14][0]*lbpar.rho_lb_units;
//      lbfluid[1][15][k] = - lbmodel.coeff[15][0]*lbpar.rho_lb_units;
//      lbfluid[1][18][k] = - lbmodel.coeff[18][0]*lbpar.rho_lb_units;

    }

  }
  
  /* top-down sweep */
//  for (k=lblattice.halo_grid_volume-lblattice.halo_offset-1;k>=0;k--) {
//
//    if (lbfields[k].boundary) {
//
//      /* bounce back to higher indices */
//      lbfluid[1][reverse[6]][k-next[6]]   = lbfluid[1][6][k];
//      lbfluid[1][reverse[12]][k-next[12]] = lbfluid[1][12][k];
//      lbfluid[1][reverse[13]][k-next[13]] = lbfluid[1][13][k];
//      lbfluid[1][reverse[16]][k-next[16]] = lbfluid[1][16][k];
//      lbfluid[1][reverse[17]][k-next[17]] = lbfluid[1][17][k];
//      
//      /* delete populations in the wall */
//      lbfluid[1][6][k]  = - lbmodel.coeff[6][0]*lbpar.rho_lb_units;
//      lbfluid[1][12][k] = - lbmodel.coeff[12][0]*lbpar.rho_lb_units;
//      lbfluid[1][13][k] = - lbmodel.coeff[13][0]*lbpar.rho_lb_units;
//      lbfluid[1][16][k] = - lbmodel.coeff[16][0]*lbpar.rho_lb_units;
//      lbfluid[1][17][k] = - lbmodel.coeff[17][0]*lbpar.rho_lb_units;
//
//    }
//
//  }
#else
#error Bounce back boundary conditions are only implemented for PUSH scheme!
#endif
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
#ifndef PULL
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
  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++) {

    if (lbfields[k].boundary) {

      /* reflect to lower indices */
      lbfluid[1][reflect[5]][k-next[5]] = lbfluid[1][5][k];
      lbfluid[1][reflect[11]][k-next[11]] = lbfluid[1][11][k];
      lbfluid[1][reflect[14]][k-next[14]] = lbfluid[1][14][k];
      lbfluid[1][reflect[15]][k-next[15]] = lbfluid[1][15][k];
      lbfluid[1][reflect[18]][k-next[18]] = lbfluid[1][18][k];

      /* delete populations in the wall */
      lbfluid[1][5][k] = -lbmodel.coeff[5][0]*lbpar.rho_lb_units;
      lbfluid[1][11][k] = -lbmodel.coeff[11][0]*lbpar.rho_lb_units;
      lbfluid[1][14][k] = -lbmodel.coeff[14][0]*lbpar.rho_lb_units;
      lbfluid[1][15][k] = -lbmodel.coeff[15][0]*lbpar.rho_lb_units;
      lbfluid[1][18][k] = -lbmodel.coeff[18][0]*lbpar.rho_lb_units;

    }
    
  }

  /* top-down sweep */
  for (k=lblattice.halo_grid_volume-lblattice.halo_offset-1;k>=0;k--) {

    if (lbfields[k].boundary) {

      /* reflect to higher indices */
      lbfluid[1][reflect[6]][k-next[6]] = lbfluid[1][6][k];
      lbfluid[1][reflect[12]][k-next[12]] = lbfluid[1][12][k];
      lbfluid[1][reflect[13]][k-next[13]] = lbfluid[1][13][k];
      lbfluid[1][reflect[16]][k-next[16]] = lbfluid[1][16][k];
      lbfluid[1][reflect[17]][k-next[17]] = lbfluid[1][17][k];

      /* delete populations in the wall */
      lbfluid[1][6][k] = -lbmodel.coeff[6][0]*lbpar.rho_lb_units;
      lbfluid[1][12][k] = -lbmodel.coeff[12][0]*lbpar.rho_lb_units;
      lbfluid[1][13][k] = -lbmodel.coeff[13][0]*lbpar.rho_lb_units;
      lbfluid[1][16][k] = -lbmodel.coeff[16][0]*lbpar.rho_lb_units;
      lbfluid[1][17][k] = -lbmodel.coeff[17][0]*lbpar.rho_lb_units;

    }
    
  }
#else
#error Specular reflections are only implemented for PUSH scheme!
#endif
#else
#error Specular reflections are only implemented for D3Q19!
#endif
}

/* slip reflection boundary condition */
MDINLINE void lb_slip_reflection() {
#if 0 //problem with slip_pref
#ifdef D3Q19
#ifndef PULL
  int k;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[19];
  next[0]  =   0;                     // ( 0, 0, 0) =
  next[1]  =   1;                     // ( 1, 0, 0) +
  next[2]  = - 1;                     // (-1, 0, 0)
  next[3]  =   yperiod;               // ( 0, 1, 0) +
  next[4]  = - yperiod;               // ( 0,-1, 0)
  next[5]  =   zperiod;               // ( 0, 0, 1) +
  next[6]  = - zperiod;               // ( 0, 0,-1)
  next[7]  =   (1+yperiod);           // ( 1, 1, 0) +
  next[8]  = - (1+yperiod);           // (-1,-1, 0)
  next[9]  =   (1-yperiod);           // ( 1,-1, 0) 
  next[10] = - (1-yperiod);           // (-1, 1, 0) +
  next[11] =   (1+zperiod);           // ( 1, 0, 1) +
  next[12] = - (1+zperiod);           // (-1, 0,-1)
  next[13] =   (1-zperiod);           // ( 1, 0,-1)
  next[14] = - (1-zperiod);           // (-1, 0, 1) +
  next[15] =   (yperiod+zperiod);     // ( 0, 1, 1) +
  next[16] = - (yperiod+zperiod);     // ( 0,-1,-1)
  next[17] =   (yperiod-zperiod);     // ( 0, 1,-1)
  next[18] = - (yperiod-zperiod);     // ( 0,-1, 1) +
  //int reflect[] = { 0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 13, 14, 11, 12, 17, 18, 15, 16 };
  //int reverse[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };

  double s = lb_boundary_par.slip_pref; //doesnt work like that any more. slip_pref has to be coded in lbfields (georg, 03.08.10)
  double r = 1.0 - s;

  double **n = lbfluid[1];

  /* bottom-up sweep */
  for (k=lblattice.halo_offset; k<lblattice.halo_grid_volume; k++) {

    if (lbfields[k].boundary) {      

      /* slip reflect to lower indices */
      n[6][k-zperiod]  = n[5][k];
      n[12][k-zperiod] = s*n[14][k] + r*n[11][k-zperiod+next[11]];
      n[13][k-zperiod] = s*n[11][k] + r*n[14][k-zperiod+next[14]];
      n[16][k-zperiod] = s*n[18][k] + r*n[15][k-zperiod+next[15]];
      n[17][k-zperiod] = s*n[15][k] + r*n[18][k-zperiod+next[18]];

      /* delete upgoing velocities (take care to delete the right ones!) */
      n[5][k]                   = -lbmodel.coeff[5][0]*lbpar.rho_lb_units;
      n[11][k]                  = -lbmodel.coeff[11][0]*lbpar.rho_lb_units;
      n[14][k-zperiod+next[14]] = -lbmodel.coeff[14][0]*lbpar.rho_lb_units;
      n[15][k]                  = -lbmodel.coeff[15][0]*lbpar.rho_lb_units;
      n[18][k-zperiod+next[18]] = -lbmodel.coeff[18][0]*lbpar.rho_lb_units;

    }

  }

  /* top-down sweep */
  for (k=lblattice.halo_grid_volume-lblattice.halo_offset-1; k>=0; k--) {

    if (lbfields[k].boundary) {

      /* slip reflect from higher indices */
      n[5][k+zperiod]  = n[6][k];
      n[11][k+zperiod] = s*n[13][k] + r*n[12][k+zperiod+next[12]];
      n[14][k+zperiod] = s*n[12][k] + r*n[13][k+zperiod+next[13]];
      n[15][k+zperiod] = s*n[17][k] + r*n[16][k+zperiod+next[16]];
      n[18][k+zperiod] = s*n[16][k] + r*n[17][k+zperiod+next[17]];

      /* delete downgoing velocities (take care to delete the right ones!) */
      n[6][k]                   = -lbmodel.coeff[6][0]*lbpar.rho_lb_units;
      n[12][k]                  = -lbmodel.coeff[12][0]*lbpar.rho_lb_units;
      n[13][k+zperiod+next[13]] = -lbmodel.coeff[13][0]*lbpar.rho_lb_units;
      n[16][k]                  = -lbmodel.coeff[16][0]*lbpar.rho_lb_units;
      n[17][k+zperiod+next[17]] = -lbmodel.coeff[17][0]*lbpar.rho_lb_units;

    }

  }
#else
#error Slip reflections are only implemented for PUSH scheme!
#endif
#else
#error Slip reflections are only implemented for D3Q19!
#endif
#endif //if 0
}

MDINLINE void lb_local_reflection(int index) {

  //int reverse[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 }; /* bounce back */
  int reverse[] = { 0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 13, 14, 11, 12, 17, 18, 15, 16 }; /* specular reflection */

  if (lbfields[index].boundary > 0) {
    lbfluid[0][5][index]  = lbfluid[0][reverse[5]][index];
    lbfluid[0][11][index] = lbfluid[0][reverse[11]][index];
    lbfluid[0][14][index] = lbfluid[0][reverse[14]][index];
    lbfluid[0][15][index] = lbfluid[0][reverse[15]][index];
    lbfluid[0][18][index] = lbfluid[0][reverse[18]][index];    
  } else {
    lbfluid[0][6][index]  = lbfluid[0][reverse[6]][index];
    lbfluid[0][12][index] = lbfluid[0][reverse[12]][index];
    lbfluid[0][13][index] = lbfluid[0][reverse[13]][index];
    lbfluid[0][16][index] = lbfluid[0][reverse[16]][index];
    lbfluid[0][17][index] = lbfluid[0][reverse[17]][index];
  }

}


MDINLINE void lb_modified_bounce_back(int index) {

  int i;
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
  //int reverse[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };
  int reverse[] = { 0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 13, 14, 11, 12, 17, 18, 15, 16 };

  double rho, j[3], n_eq[19], avg_rho=lbpar.rho_lb_units;
  double (*coeff)[4]=lbmodel.coeff, (*c)[3]=lbmodel.c;

  lb_calc_local_j(index,j);
  lb_calc_local_rho(index,&rho);

  if (lbfields[index].boundary > 0) {
//  fprintf(stderr,"%+d rho=%f j=(%f,%f,%f)\n",lbfields[index].boundary,rho,j[0],j[1],j[2]);
  }

  ///* precollision unknown distributions */
  //if (lbfields[index].boundary > 0) {
  //  lbfluid[0][5][index]  = lbfluid[0][reverse[5]][index];
  //  lbfluid[0][11][index] = lbfluid[0][reverse[11]][index];
  //  lbfluid[0][14][index] = lbfluid[0][reverse[14]][index];
  //  lbfluid[0][15][index] = lbfluid[0][reverse[15]][index];
  //  lbfluid[0][18][index] = lbfluid[0][reverse[18]][index];    
  //} else {
  //  lbfluid[0][6][index]  = lbfluid[0][reverse[6]][index];
  //  lbfluid[0][12][index] = lbfluid[0][reverse[12]][index];
  //  lbfluid[0][13][index] = lbfluid[0][reverse[13]][index];
  //  lbfluid[0][16][index] = lbfluid[0][reverse[16]][index];
  //  lbfluid[0][17][index] = lbfluid[0][reverse[17]][index];
  //}

  if (lbfields[index].boundary > 0) {
    fprintf(stderr,"pre coll n=(");
    for (i=0;i<lbmodel.n_veloc;i++) {
      fprintf(stderr," %.9f",lbfluid[0][i][index]+coeff[i][0]);
    }
    fprintf(stderr," \n");
  }

  /* lb_calc_local_j addiert schon die Kraft dazu!!!! */

  lb_calc_local_j(index,j);
  lb_calc_local_rho(index,&rho);

  if (lbfields[index].boundary > 0) {
//  fprintf(stderr,"%+d rho=%f j=(%f,%f,%f)\n",lbfields[index].boundary,rho,j[0],j[1],j[2]);
  }

  double u[3], f[3], C[6];

  /* collisions */
  lb_calc_local_rho(index, &rho);
  lb_calc_local_j(index, j);

  f[0] = lbpar.ext_force[0];
  f[1] = lbpar.ext_force[1];
  f[2] = lbpar.ext_force[2];

  //j[0] += 0.5*f[0];
  //j[1] += 0.5*f[1];
  //j[2] += 0.5*f[2];

  u[0] = (j[0]+0.5*f[0])/rho;
  u[1] = (j[1]+0.5*f[1])/rho;
  u[2] = (j[2]+0.5*f[2])/rho;

  //fprintf(stderr,"rho=%f j=(%f,%f,%f)\n",rho,j[0],j[1],j[2]);

  for (i=0;i<lbmodel.n_veloc;i++) {
    n_eq[i] =  coeff[i][0] * (rho-avg_rho);
    n_eq[i] += coeff[i][1] * scalar(j,c[i]);
    n_eq[i] += coeff[i][2] * SQR(scalar(u,c[i])) * rho;
    n_eq[i] += coeff[i][3] * scalar(u,u) * rho;
  }

  double gamma_shear = 1. - 2./(6.*lbpar.viscosity*lbpar.tau/SQR(lbpar.agrid)+1.);

  for (i=0; i<lbmodel.n_veloc; i++) {
    lbfluid[0][i][index] += (gamma_shear-1.)*(lbfluid[0][i][index]-n_eq[i]);
  }

  lb_calc_local_j(index,j);
  lb_calc_local_rho(index,&rho);

  if (lbfields[index].boundary > 0) {
//  fprintf(stderr,"%+d rho=%f j=(%f,%f,%f) f=(%f,%f,%f)\n",lbfields[index].boundary,rho,j[0],j[1],j[2],f[0],f[1],f[2]);
  }

  /* forcing */
  C[0] = (1.+gamma_shear)*u[0]*f[0];// + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[2] = (1.+gamma_shear)*u[1]*f[1];// + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[5] = (1.+gamma_shear)*u[2]*f[2];// + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[1] = 1./2.*(1.+gamma_shear)*(u[0]*f[1]+u[1]*f[0]);
  C[3] = 1./2.*(1.+gamma_shear)*(u[0]*f[2]+u[2]*f[0]);
  C[4] = 1./2.*(1.+gamma_shear)*(u[1]*f[2]+u[2]*f[1]);

  double tmp =   C[0]*c[i][0]*c[i][0] 
               + (2.*C[1]*c[i][0]+C[2]*c[i][1])*c[i][1]
               + (2.*(C[3]*c[i][0]+C[4]*c[i][1])+C[5]*c[i][2])*c[i][2];
  double trace = C[0]+C[2]+C[3];

  for (i=0; i<lbmodel.n_veloc; i++) {
    lbfluid[0][i][index] += coeff[i][1]*scalar(f,c[i]);
    lbfluid[0][i][index] += coeff[i][2]*tmp;
    lbfluid[0][i][index] += coeff[i][3]*trace;
  }  

  //f[0] = - lb_boundary_par.slip_pref * j[0]/rho;
  //f[1] = - lb_boundary_par.slip_pref * j[1]/rho;
  //f[2] = - lb_boundary_par.slip_pref * j[2]/rho;

  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  lbfluid[0][i][index] += coeff[i][1]*scalar(f,c[i]);
  //}  

  //if (lbfields[index].boundary > 0) {
  //  lbfluid[0][11][index] +=  0.5 * f[0];
  //  lbfluid[0][14][index] += -0.5 * f[0];
  //  lbfluid[0][15][index] +=  0.5 * f[1];
  //  lbfluid[0][18][index] += -0.5 * f[1];
  //} else {
  //  lbfluid[0][12][index] += -0.5 * f[0];
  //  lbfluid[0][13][index] +=  0.5 * f[0];
  //  lbfluid[0][16][index] += -0.5 * f[1];
  //  lbfluid[0][17][index] +=  0.5 * f[1];
  //}

  lb_calc_local_j(index,j);
  lb_calc_local_rho(index,&rho);

  if (lbfields[index].boundary > 0) {
//  fprintf(stderr,"%+d rho=%f j=(%f,%f,%f) f=(%f,%f,%f)\n",lbfields[index].boundary,rho,j[0],j[1],j[2],f[0],f[1],f[2]);
  }

  if (lbfields[index].boundary > 0) {
  fprintf(stderr,"pre push n=(");
  for (i=0;i<lbmodel.n_veloc;i++) {
    fprintf(stderr," %.9f",lbfluid[0][i][index]+coeff[i][0]);
  }
  fprintf(stderr," \n");
  }

  /* push */
  lbfluid[1][0][index+next[0]]  = lbfluid[0][0][index];
  lbfluid[1][1][index+next[1]]  = lbfluid[0][1][index];
  lbfluid[1][2][index+next[2]]  = lbfluid[0][2][index];
  lbfluid[1][3][index+next[3]]  = lbfluid[0][3][index];
  lbfluid[1][4][index+next[4]]  = lbfluid[0][4][index];
  lbfluid[1][7][index+next[7]]  = lbfluid[0][7][index];
  lbfluid[1][8][index+next[8]]  = lbfluid[0][8][index];
  lbfluid[1][9][index+next[9]]  = lbfluid[0][9][index];
  lbfluid[1][10][index+next[10]] = lbfluid[0][10][index];
  
  if (lbfields[index].boundary > 0) {
    lbfluid[1][5][index+next[5]]   = lbfluid[0][5][index];
    lbfluid[1][11][index+next[11]] = lbfluid[0][11][index];
    lbfluid[1][14][index+next[14]] = lbfluid[0][14][index];
    lbfluid[1][15][index+next[15]] = lbfluid[0][15][index];
    lbfluid[1][18][index+next[18]] = lbfluid[0][18][index];

    //lbfluid[1][6][index]  = lbfluid[0][reverse[6]][index-next[6]];
    //lbfluid[1][12][index] = lbfluid[0][reverse[12]][index-next[12]];
    //lbfluid[1][13][index] = lbfluid[0][reverse[13]][index-next[13]];
    //lbfluid[1][16][index] = lbfluid[0][reverse[16]][index-next[16]];
    //lbfluid[1][17][index] = lbfluid[0][reverse[17]][index-next[17]];
  } else {
    lbfluid[1][6][index+next[6]]   = lbfluid[0][6][index];
    lbfluid[1][12][index+next[12]] = lbfluid[0][12][index];
    lbfluid[1][13][index+next[13]] = lbfluid[0][13][index];
    lbfluid[1][16][index+next[16]] = lbfluid[0][16][index];
    lbfluid[1][17][index+next[17]] = lbfluid[0][17][index];

    //lbfluid[1][5][index]  = lbfluid[0][reverse[5]][index-next[5]];
    //lbfluid[1][11][index] = lbfluid[0][reverse[11]][index-next[11]];
    //lbfluid[1][14][index] = lbfluid[0][reverse[14]][index-next[14]];
    //lbfluid[1][15][index] = lbfluid[0][reverse[15]][index-next[15]];
    //lbfluid[1][18][index] = lbfluid[0][reverse[18]][index-next[18]];
  }

  //if (lbfields[index].boundary > 0) {
  //fprintf(stderr,"post push n=(");
  //for (i=0;i<lbmodel.n_veloc;i++) {
  //  fprintf(stderr," %.9f",lbfluid[1][i][index]+coeff[i][0]);
  //}
  //fprintf(stderr," \n");
  //}

}

MDINLINE void lb_boundary_push(index_t index) {

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  index_t next[19];
  next[0]  = index;
  next[1]  = index + 1;
  next[2]  = index - 1;
  next[3]  = index + yperiod;
  next[4]  = index - yperiod;
  next[5]  = index + zperiod;
  next[6]  = index - zperiod;
  next[7]  = index + (1 + yperiod);
  next[8]  = index - (1 + yperiod);
  next[9]  = index + (1 - yperiod);
  next[10] = index - (1 - yperiod);
  next[11] = index + (1 + zperiod);
  next[12] = index - (1 + zperiod);
  next[13] = index + (1 - zperiod);
  next[14] = index - (1 - zperiod);
  next[15] = index + (yperiod + zperiod);
  next[16] = index - (yperiod + zperiod);
  next[17] = index + (yperiod - zperiod);
  next[18] = index - (yperiod - zperiod);
  
  lbfluid[1][0][next[0]] = lbfluid[0][0][index];
  lbfluid[1][1][next[1]] = lbfluid[0][1][index];
  lbfluid[1][2][next[2]] = lbfluid[0][2][index];
  lbfluid[1][3][next[3]] = lbfluid[0][3][index];
  lbfluid[1][4][next[4]] = lbfluid[0][4][index];
  lbfluid[1][5][next[5]] = lbfluid[0][5][index];
  lbfluid[1][6][next[6]] = lbfluid[0][6][index];
  lbfluid[1][7][next[7]] = lbfluid[0][7][index];
  lbfluid[1][8][next[8]] = lbfluid[0][8][index];
  lbfluid[1][9][next[9]] = lbfluid[0][9][index];
  lbfluid[1][10][next[10]] = lbfluid[0][10][index];
  lbfluid[1][11][next[11]] = lbfluid[0][11][index];
  lbfluid[1][12][next[12]] = lbfluid[0][12][index];
  lbfluid[1][13][next[13]] = lbfluid[0][13][index];
  lbfluid[1][14][next[14]] = lbfluid[0][14][index];
  lbfluid[1][15][next[15]] = lbfluid[0][15][index];
  lbfluid[1][16][next[16]] = lbfluid[0][16][index];
  lbfluid[1][17][next[17]] = lbfluid[0][17][index];
  lbfluid[1][18][next[18]] = lbfluid[0][18][index];

}

MDINLINE void lb_boundary_forces(index_t index, double *mode) {
#if 0 //problem with slip_pref (georg, 03.08.10)
  int i;
  double (*coeff)[4]=lbmodel.coeff, (*c)[3]=lbmodel.c;
  double rho=mode[0]+lbpar.rho_lb_units, *j=&mode[1], *f=lbfields[index].force;
  double u[3], C[6];

  lb_calc_local_j(index,j);

  f[0] = lbpar.ext_force[0];
  f[1] = lbpar.ext_force[1];
  f[2] = lbpar.ext_force[2];

//  fprintf(stderr,"rho=%f j=(%f,%f,%f) f=(%f,%f,%f)\n",rho,j[0],j[1],j[2],f[0],f[1],f[2]);

  f[0] += -lb_boundary_par.slip_pref * (j[0]+0.5*f[0])/rho;
  f[1] += -lb_boundary_par.slip_pref * (j[1]+0.5*f[1])/rho;
  f[2] += -lb_boundary_par.slip_pref * (j[2]+0.5*f[2])/rho;

  u[0] = (j[0]+0.5*f[0])/rho;
  u[1] = (j[1]+0.5*f[1])/rho;
  u[2] = (j[2]+0.5*f[1])/rho;

  //lb_calc_local_j(index,j);
  //fprintf(stderr,"rho=%f j=(%f,%f,%f) f=(%f,%f,%f)\n",rho,j[0],j[1],j[2],f[0],f[1],f[2]);

  //if (lbfields[index].boundary < 0) {
  //
  //  lbfluid[0][11][index] +=  0.5 * f[0];
  //  lbfluid[0][14][index] += -0.5 * f[0];
  //  lbfluid[0][15][index] +=  0.5 * f[1];
  //  lbfluid[0][18][index] += -0.5 * f[1];    
  //
  //} else {
  //
  //  lbfluid[0][12][index] += -0.5 * f[0];
  //  lbfluid[0][13][index] +=  0.5 * f[0];
  //  lbfluid[0][16][index] += -0.5 * f[1];
  //  lbfluid[0][17][index] +=  0.5 * f[1];
  //  
  //}

  double gamma_shear = 1. - 2./(6.*lbpar.viscosity*lbpar.tau/SQR(lbpar.agrid)+1.);
  double gamma_bulk = 0.0;

  /* forcing */
  C[0] = (1.+gamma_shear)*u[0]*f[0] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[2] = (1.+gamma_shear)*u[1]*f[1] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[5] = (1.+gamma_shear)*u[2]*f[2] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  C[1] = 1./2.*(1.+gamma_shear)*(u[0]*f[1]+u[1]*f[0]);
  C[3] = 1./2.*(1.+gamma_shear)*(u[0]*f[2]+u[2]*f[0]);
  C[4] = 1./2.*(1.+gamma_shear)*(u[1]*f[2]+u[2]*f[1]);

  double tmp =   C[0]*c[i][0]*c[i][0] 
               + (2.*C[1]*c[i][0]+C[2]*c[i][1])*c[i][1]
               + (2.*(C[3]*c[i][0]+C[4]*c[i][1])+C[5]*c[i][2])*c[i][2];
  double trace = C[0]+C[2]+C[3];

  for (i=0; i<lbmodel.n_veloc; i++) {
    lbfluid[0][i][index] += coeff[i][1]*scalar(f,c[i]);
    lbfluid[0][i][index] += coeff[i][2]*tmp;
    lbfluid[0][i][index] += coeff[i][3]*trace;
  }  

  //if (lbfields[index].boundary > 0) {
  //  lbfluid[0][11][index] +=  0.5 * f[0];
  //  lbfluid[0][14][index] += -0.5 * f[0];
  //  lbfluid[0][15][index] +=  0.5 * f[1];
  //  lbfluid[0][18][index] += -0.5 * f[1];
  //} else {
  //  lbfluid[0][12][index] += -0.5 * f[0];
  //  lbfluid[0][13][index] +=  0.5 * f[0];
  //  lbfluid[0][16][index] += -0.5 * f[1];
  //  lbfluid[0][17][index] +=  0.5 * f[1];
  //}

  //{
  //  lbfluid[0][1][index] +=  0.5 * f[0];
  //  lbfluid[0][2][index] += -0.5 * f[0];
  //  lbfluid[0][3][index] +=  0.5 * f[1];
  //  lbfluid[0][4][index] += -0.5 * f[1];
  //}

  lb_calc_local_j(index,j);
  
//  fprintf(stderr,"rho=%f j=(%f,%f,%f) f=(%f,%f,%f)\n",rho,j[0],j[1],j[2],f[0],f[1],f[2]);
#endif //if 0
}

MDINLINE void lb_boundary_equilibrium(int index, double rho, double *v, double *pi, int df) {
  int i;

#ifdef D3Q19
  const double A = 1./6.;
  const double P = 5./36.;
  const double Z = 5./108.; // 11./108.; war falsch!

  const double w[19] = { 7./18.,
                         1./12., 1./12., 1./12., 1./12., 1./18., 1./18.,
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36. };

  double (*c)[3] = lbmodel.c;
      
  double chi1, chi2, lambda1[3], lambda2[3], lsum[3];
  double c0;

  lambda1[0] = v[0] / lbmodel.c_sound_sq;
  lambda1[1] = v[1] / lbmodel.c_sound_sq;
  lambda1[2] = (v[2])/P;
  
  chi1 = -A*lambda1[2];

  lambda2[0] = 0.0;
  lambda2[1] = 0.0;
  lambda2[2] = -Z/P*SQR(lambda1[2]);

  chi2 = -A*lambda2[2] - 0.5*P*SQR(lambda1[2]) - 1./6.*(SQR(lambda1[0])+SQR(lambda1[1]));

  c0 = chi1;// + 0.5*SQR(chi1) + chi2;

  lsum[0] = lambda1[0];// + lambda2[0] + lambda1[0]*chi1;
  lsum[1] = lambda1[1];// + lambda2[1] + lambda1[1]*chi1;
  lsum[2] = lambda1[2];// + lambda2[2] + lambda1[2]*chi1;

  for (i=0; i<lbmodel.n_veloc; i++) {

    if (df*scalar(c[i],lbfields[index].nvec) >= 0.0) {
      
      lbfluid[0][i][index] =   w[i] * rho
			     + w[i] * rho * ( c0 
					      + scalar(lsum,c[i]) 
					      + 0.5*SQR(scalar(lambda1,c[i])) )
	- lbmodel.coeff[i][0]*lbpar.rho_lb_units;

    } else {

      lbfluid[0][i][index] = -lbmodel.coeff[i][0]*lbpar.rho_lb_units;

    }

  }

  if(index == get_linear_index(1,1,1,lblattice.halo_grid)) {
    double rho=0.0, j[3], pi[6];
    for (i=0;i<lbmodel.n_veloc;i++) {
	rho += lbfluid[0][i][index];
	j[0] += lbfluid[0][i][index]*c[i][0];
	j[1] += lbfluid[0][i][index]*c[i][1];
	j[2] += lbfluid[0][i][index]*c[i][2];
	pi[0] += lbfluid[0][i][index]*c[i][0]*c[i][0];
	pi[1] += lbfluid[0][i][index]*c[i][0]*c[i][1];
	pi[2] += lbfluid[0][i][index]*c[i][1]*c[i][1];
	pi[3] += lbfluid[0][i][index]*c[i][0]*c[i][2];
	pi[4] += lbfluid[0][i][index]*c[i][1]*c[i][2];
	pi[5] += lbfluid[0][i][index]*c[i][2]*c[i][2];
    }   
    //double *n = lbfields[index].nvec;
    //fprintf(stderr,"nvec = (%e,%e,%e)\n",n[0],n[1],n[2]);
    //fprintf(stderr,"A = %e P = %e Z = %e\n",A,P,Z);
    //fprintf(stderr,"chi1 = %e chi2 = %f lambda1 = (%e,%e,%e) lambda2 = (%e,%e,%e)\n",chi1,chi2,lambda1[0],lambda1[1],lambda1[2],lambda2[0],lambda2[1],lambda2[2]);
    //fprintf(stderr,"rho = %e j = (%e,%e,%e) pi = (%e,%e,%e,%e,%e,%e)\n",rho,j[0],j[1],j[2],pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);
    //for (i=0; i<lbmodel.n_veloc; i++) {      
    //  fprintf(stderr,"[%d] %e\n",i,lbfluid[0][i][index]+lbmodel.coeff[i][0]*lbpar.rho_lb_units);
    //}
  }
#else
#error lb_boundary_equilibrium_precoll() is only implemented for D3Q19
#endif
}

MDINLINE void lb_boundary_calc_modes(int index, double *mode, double *pi) {

  int i;

  fprintf(stderr,"lb_boundary_calc_modes(%f)\n",lbfields[index].nvec[2]);

  //double pi[6];
  lb_calc_local_pi(index, pi);
  fprintf(stderr,"pi=(%f,%f,%f,%e,%f,%f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5],pi[6]);


  for (i=0;i<lbmodel.n_veloc;i++) {
    lbfluid[0][i][index] += lbmodel.coeff[i][0];
  }

  //fprintf(stderr,"( ");
  //for (i=0;i<lbmodel.n_veloc;i++) {
  //  fprintf(stderr,"%f ",lbfluid[0][i][index]);
  //}
  //fprintf(stderr,") %f\n",lbfields[index].nvec[2]);
  
  //lb_calc_modes(index, mode);
  
  //fprintf(stderr,"(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f)\n",mode[0],mode[1],mode[2],mode[3],mode[4],mode[5],mode[6],mode[7],mode[8],mode[9]);

  /* first we do the specular reflections */
  //if (lbfields[index].nvec[2] < 0.0) {
  //  lbfluid[0][6][index] = lbfluid[0][5][index];
  //  lbfluid[0][12][index] = lbfluid[0][14][index];
  //  lbfluid[0][13][index] = lbfluid[0][11][index];
  //  lbfluid[0][16][index] = lbfluid[0][18][index];
  //  lbfluid[0][17][index] = lbfluid[0][15][index];
  //  lbfluid[0][5][index] = 0.0;
  //  lbfluid[0][11][index] = 0.0;
  //  lbfluid[0][14][index] = 0.0;
  //  lbfluid[0][15][index] = 0.0;
  //  lbfluid[0][18][index] = 0.0;
  //} else {
  //  lbfluid[0][5][index] = lbfluid[0][6][index];
  //  lbfluid[0][11][index] = lbfluid[0][13][index];
  //  lbfluid[0][14][index] = lbfluid[0][12][index];
  //  lbfluid[0][15][index] = lbfluid[0][17][index];
  //  lbfluid[0][18][index] = lbfluid[0][16][index];
  //  lbfluid[0][6][index] = 0.0;
  //  lbfluid[0][12][index] = 0.0;
  //  lbfluid[0][13][index] = 0.0;
  //  lbfluid[0][16][index] = 0.0;
  //  lbfluid[0][17][index] = 0.0;
  //}

  //lb_calc_modes(index, mode);

  //fprintf(stderr,"(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f)\n",mode[0],mode[1],mode[2],mode[3],mode[4],mode[5],mode[6],mode[7],mode[8],mode[9]);

  double n0,n1p,n1m,n2p,n2m,n3,n4p,n4m,n5p,n5m,n6p,n6m,n7p,n7m;

  n0  = lbfluid[0][0][index];
  n1p = lbfluid[0][1][index] + lbfluid[0][2][index];
  n1m = lbfluid[0][1][index] - lbfluid[0][2][index];
  n2p = lbfluid[0][3][index] + lbfluid[0][4][index];
  n2m = lbfluid[0][3][index] - lbfluid[0][4][index];
  n4p = lbfluid[0][7][index] + lbfluid[0][8][index];
  n4m = lbfluid[0][7][index] - lbfluid[0][8][index];
  n5p = lbfluid[0][9][index] + lbfluid[0][10][index];
  n5m = lbfluid[0][9][index] - lbfluid[0][10][index];

  if (lbfields[index].nvec[2] > 0.0) {

    n3  = lbfluid[0][6][index];
    n6p = lbfluid[0][13][index] + lbfluid[0][12][index];
    n6m = lbfluid[0][13][index] - lbfluid[0][12][index];
    n7p = lbfluid[0][17][index] + lbfluid[0][16][index];
    n7m = lbfluid[0][17][index] - lbfluid[0][16][index];

  } else {

    n3  = lbfluid[0][5][index];
    n6p = lbfluid[0][11][index] + lbfluid[0][14][index];
    n6m = lbfluid[0][11][index] - lbfluid[0][14][index];
    n7p = lbfluid[0][15][index] + lbfluid[0][18][index];
    n7m = lbfluid[0][15][index] - lbfluid[0][18][index];

  }
    
  mode[0] = n0 + n1p + n2p + n3 + n4p + n5p + n6p + n7p;

  mode[1] = n1m + n4m + n5m + n6m;
  mode[2] = n2m + n4m - n5m + n7m;
  mode[3] = -lbfields[index].nvec[2]*(-n0 - n1p - n2p + 5.*n3 - n4p - n5p + 5.*(n6p+n7p))/6.;

  mode[1] = 0.0;
  mode[2] = 0.0;
  mode[3] = 0.0;

//  fprintf(stderr, "rho=%f j=(%e,%f,%f)\n",mode[0],mode[1],mode[2],mode[3]);

  const double A = lbfields[index].nvec[2] * 1./6.;
  const double P = 5./36.;
  const double Z = 5./108.; // 11./108.; war falsch!

  const double w[19] = { 7./18.,
                         1./12., 1./12., 1./12., 1./12., 1./18., 1./18.,
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36. };

  double (*c)[3] = lbmodel.c;
      
  double chi1, chi2, lambda1[3], lambda2[3], lsum[3];
  double c0;

  lambda1[0] = (mode[1]/mode[0])/lbmodel.c_sound_sq;
  lambda1[1] = (mode[2]/mode[0])/lbmodel.c_sound_sq;
  lambda1[2] = (mode[3]/mode[0])/P;
  
  chi1 = -A*lambda1[2];

  lambda2[0] = 0.0;
  lambda2[1] = 0.0;
  lambda2[2] = -Z/P*SQR(lambda1[2]);

  chi2 = -A*lambda2[2] - 0.5*P*SQR(lambda1[2]) - 1./6.*(SQR(lambda1[0])+SQR(lambda1[1]));

  lsum[0] = lambda1[0] + lambda2[0] + lambda1[0]*chi1;
  lsum[1] = lambda1[1] + lambda2[1] + lambda1[1]*chi1;
  lsum[2] = lambda1[2] + lambda2[2] + lambda1[2]*chi1;

  double n_eq[19], n_neq[19];

  //if (lbfields[index].nvec[2] > 0.0) {
  //
  //  n_eq[6]  = w[6]  * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[6]) 
  //				   + 0.5*SQR(scalar(lambda1,c[6])) );
  //  n_eq[12] = w[12] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[12]) 
  //				   + 0.5*SQR(scalar(lambda1,c[12])) );
  //  n_eq[13] = w[13] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[13]) 
  //				   + 0.5*SQR(scalar(lambda1,c[13])) );
  //  n_eq[16] = w[16] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[16]) 
  //				   + 0.5*SQR(scalar(lambda1,c[16])) );
  //  n_eq[17] = w[17] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[17]) 
  //				   + 0.5*SQR(scalar(lambda1,c[17])) );
  //
  //  n_neq[6]  = lbfluid[0][6][index] - n_eq[6];
  //  n_neq[12] = lbfluid[0][12][index] - n_eq[12];
  //  n_neq[13] = lbfluid[0][13][index] - n_eq[13];
  //  n_neq[16] = lbfluid[0][16][index] - n_eq[16];
  //  n_neq[17] = lbfluid[0][17][index] - n_eq[17];
  //
  //} else {
  //
  //  n_eq[5]  = w[5]  * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[5]) 
  //				   + 0.5*SQR(scalar(lambda1,c[5])) );
  //  n_eq[11] = w[11] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[11]) 
  //				   + 0.5*SQR(scalar(lambda1,c[11])) );
  //  n_eq[14] = w[14] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[14]) 
  //				   + 0.5*SQR(scalar(lambda1,c[14])) );
  //  n_eq[15] = w[15] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[15]) 
  //				   + 0.5*SQR(scalar(lambda1,c[15])) );
  //  n_eq[18] = w[18] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  //				   + scalar(lsum,c[18]) 
  //				   + 0.5*SQR(scalar(lambda1,c[18])) );
  //
  //  n_neq[5]  = lbfluid[0][5][index] - n_eq[5];
  //  n_neq[11] = lbfluid[0][11][index] - n_eq[11];
  //  n_neq[14] = lbfluid[0][14][index] - n_eq[14];
  //  n_neq[15] = lbfluid[0][15][index] - n_eq[15];
  //  n_neq[18] = lbfluid[0][18][index] - n_eq[18];
  //
  //}

  for (i=0; i<lbmodel.n_veloc; i++) {
    n_eq[i] = w[i] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  				   + scalar(lsum,c[i]) 
  				   + 0.5*SQR(scalar(lambda1,c[i])) );
  }
  
  for (i=0; i<lbmodel.n_veloc; i++) {
    n_neq[i] = lbfluid[0][i][index] - n_eq[i];
  }
 
  fprintf(stderr,"n_eq=(\n");
  for (i=0; i<lbmodel.n_veloc; i++) {
    fprintf(stderr, "%d %f %f %f\n",i,lbfluid[0][i][index],n_eq[i],n_neq[i]);
  }
  fprintf(stderr,")\n");

  int reverse[19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };

  if (lbfields[index].nvec[2] > 0.0) {
 
    pi[0] =   lbfluid[0][1][index]  + lbfluid[0][2][index]  
            + lbfluid[0][7][index]  + lbfluid[0][8][index]  
            + lbfluid[0][9][index]  + lbfluid[0][10][index]
            + lbfluid[0][11][index] + lbfluid[0][12][index] 
            + lbfluid[0][13][index] + lbfluid[0][14][index];
    pi[2] =   lbfluid[0][3][index]  + lbfluid[0][4][index]
            + lbfluid[0][7][index]  + lbfluid[0][8][index]  
            + lbfluid[0][9][index]  + lbfluid[0][10][index]
            + lbfluid[0][15][index] + lbfluid[0][16][index] 
            + lbfluid[0][17][index] + lbfluid[0][18][index];
    pi[1] =   lbfluid[0][7][index]  - lbfluid[0][9][index] 
            + lbfluid[0][8][index]  - lbfluid[0][10][index];
    pi[5] =   n_eq[5]  + n_neq[reverse[5]]  + lbfluid[0][6][index]
            + n_eq[11] + n_neq[reverse[11]] + lbfluid[0][12][index] 
            + n_eq[14] + n_neq[reverse[14]] + lbfluid[0][13][index]
            + n_eq[15] + n_neq[reverse[15]] + lbfluid[0][16][index] 
            + n_eq[18] + n_neq[reverse[18]] + lbfluid[0][17][index];
    pi[3] =   n_eq[11] + n_neq[reverse[11]] + lbfluid[0][12][index] 
            - n_eq[14] - n_neq[reverse[14]] - lbfluid[0][13][index];
    pi[4] =   n_eq[15] + n_neq[reverse[15]] + lbfluid[0][16][index]
            - n_eq[18] - n_neq[reverse[18]] - lbfluid[0][17][index];

  } else {

    pi[0] =   lbfluid[0][1][index]  + lbfluid[0][2][index]  
            + lbfluid[0][7][index]  + lbfluid[0][8][index]  
            + lbfluid[0][9][index]  + lbfluid[0][10][index]
            + lbfluid[0][11][index] + lbfluid[0][12][index] 
            + lbfluid[0][13][index] + lbfluid[0][14][index];
    pi[2] =   lbfluid[0][3][index]  + lbfluid[0][4][index]
            + lbfluid[0][7][index]  + lbfluid[0][8][index]  
            + lbfluid[0][9][index]  + lbfluid[0][10][index]
            + lbfluid[0][15][index] + lbfluid[0][16][index] 
            + lbfluid[0][17][index] + lbfluid[0][18][index];
    pi[1] =   lbfluid[0][7][index]  - lbfluid[0][9][index] 
            + lbfluid[0][8][index]  - lbfluid[0][10][index];
    pi[5] =   n_eq[6]  + n_neq[reverse[6]]  + lbfluid[0][5][index]
            + n_eq[12] + n_neq[reverse[12]] + lbfluid[0][11][index] 
            + n_eq[13] + n_neq[reverse[13]] + lbfluid[0][14][index]
            + n_eq[16] + n_neq[reverse[16]] + lbfluid[0][15][index] 
            + n_eq[17] + n_neq[reverse[17]] + lbfluid[0][18][index];
    pi[3] =   n_eq[12] + n_neq[reverse[12]] + lbfluid[0][11][index] 
            - n_eq[13] - n_neq[reverse[13]] - lbfluid[0][14][index];
    pi[4] =   n_eq[16] + n_neq[reverse[16]] + lbfluid[0][15][index]
            - n_eq[17] - n_neq[reverse[17]] - lbfluid[0][18][index];

  }

  fprintf(stderr,"pi=(%f,%f,%f,%e,%f,%f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);

  fprintf(stderr,"m7_pre = %e\n",lbfluid[0][11][index]-lbfluid[0][14][index]);
  fprintf(stderr,"pi_^neq_xz = %e\n",2.*(n_neq[11]-n_neq[14]));
  fprintf(stderr,"pi_xz = %e\n",n_eq[12]+n_neq[reverse[12]]-n_eq[13]-n_neq[reverse[13]]);
  fprintf(stderr,"pi_xz - 1./6.*j_x = %e\n",n_eq[12]+n_neq[reverse[12]]-n_eq[13]-n_neq[reverse[13]]-1./6.*mode[1]);

  mode[4] = (-2.*n0 + n1p + n2p - 2.*n3 + 4.*(n4p+n5p) + n6p + n7p)/3.;
  mode[5] = n1p - n2p + n6p - n7p;
  mode[6] = n4p - n5p;

  mode[9]  = (-2.*n1m + 3.*(n4m+n5m))/5.;
  mode[10] = (-2.*n2m + 3.*(n4m-n5m))/5.;
  mode[13] = (12.*n0 - 28.*(n1p+n2p) + 42.*(n4p+n5p))/55.;

  if (lbfields[index].nvec[2] > 0.0) {

    n3 = n_eq[5] + n_neq[reverse[5]];
    n6p = n_eq[11] + n_eq[14] + n_neq[reverse[11]] + n_neq[reverse[14]];
    n6m = n_eq[11] - n_eq[14] + n_neq[reverse[11]] - n_neq[reverse[14]];
    n7p = n_eq[15] + n_eq[18] + n_neq[reverse[15]] + n_neq[reverse[18]];
    n7m = n_eq[15] - n_eq[18] + n_neq[reverse[15]] - n_neq[reverse[18]];

    mode[3] = -(-n0 - n1p - n2p + 5.*n3 - n4p - n5p + 5.*(n6p+n7p))/6.;

    //mode[7] = -(-n1m - n4m - n5m + 5.*n6m)/6.;
    //mode[8] = -(-n2m - n4m + n5m + 5.*n7m)/6.;

    mode[7] = n6m + A * mode[1];
    mode[8] = n7m + A * mode[2];

    mode[11] = -(2.*n0 - n1p - n2p - 22.*n3 - 4.*(n4p+n5p) + 11.*(n6p+n7p))/36.;
    mode[12] = -(-n1p + n2p +3.*(n6p-n7p))/4.;

  } else {

    n3 = n_eq[6] + n_neq[reverse[6]];
    n6p = n_eq[13] + n_eq[12] + n_neq[reverse[13]] + n_neq[reverse[12]];
    n6m = n_eq[13] - n_eq[12] + n_neq[reverse[13]] - n_neq[reverse[12]];
    n7p = n_eq[17] + n_eq[16] + n_neq[reverse[17]] + n_neq[reverse[16]];
    n7m = n_eq[17] - n_eq[16] + n_neq[reverse[17]] - n_neq[reverse[16]];

    mode[3] = (-n0 - n1p - n2p + 5.*n3 - n4p - n5p + 5.*(n6p+n7p))/6.;

    //mode[7] = (-n1m - n4m - n5m + 5.*n6m)/6.;
    //mode[8] = (-n2m - n4m + n5m + 5.*n7m)/6.;

    mode[7] = - n6m + A*mode[1];
    mode[8] = - n7m + A*mode[2];

    mode[11] = (2.*n0 - n1p - n2p - 22.*n3 - 4.*(n4p+n5p) + 11.*(n6p+n7p))/36.;
    mode[12] = (-n1p + n2p +3.*(n6p-n7p))/4.;

  }

  fprintf(stderr,"modes=(%f,%e,%f,%f,%e,%e,%e,%e,%e,%f,%f,%f,%f,%f) %f\n",mode[0],mode[1],mode[2],mode[3],mode[4],mode[5],mode[6],mode[7],mode[8],mode[9],mode[10],mode[11],mode[12],mode[13],lbfields[index].nvec[2]);

}

MDINLINE void lb_boundary_bb_neq_BGK(index_t index, double *mode) {
#if 0 //problem with slip pref (georg, 03.08.10)
  int i;

//  fprintf(stderr,"lb_boundary_bb_neq_BGK(%+.0f)\n",lbfields[index].nvec[2]);

  lb_calc_local_fields(index, &mode[0], &mode[1], NULL);

//  fprintf(stderr, "rho=%f j=(%e,%f,%f)\n",mode[0],mode[1],mode[2],mode[3]);

  mode[1] = lb_boundary_par.slip_pref;
  mode[2] = 0.0;
  mode[3] = 0.0;

  for (i=0;i<lbmodel.n_veloc;i++) {
    lbfluid[0][i][index] += lbmodel.coeff[i][0];
  }

  const double A = lbfields[index].nvec[2] * 1./6.;
  const double P = 5./36.;
  const double Z = 5./108.; // 11./108.; war falsch!

  const double w[19] = { 7./18.,
                         1./12., 1./12., 1./12., 1./12., 1./18., 1./18.,
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36. };

  double (*c)[3] = lbmodel.c;
      
  double chi1, chi2, lambda1[3], lambda2[3], lsum[3];
  double c0;

  lambda1[0] = (mode[1]/mode[0])/lbmodel.c_sound_sq;
  lambda1[1] = (mode[2]/mode[0])/lbmodel.c_sound_sq;
  lambda1[2] = (mode[3]/mode[0])/P;
  
  chi1 = -A*lambda1[2];

  lambda2[0] = 0.0;
  lambda2[1] = 0.0;
  lambda2[2] = -Z/P*SQR(lambda1[2]);

  chi2 = -A*lambda2[2] - 0.5*P*SQR(lambda1[2]) - 1./6.*(SQR(lambda1[0])+SQR(lambda1[1]));

  lsum[0] = lambda1[0] + lambda2[0] + lambda1[0]*chi1;
  lsum[1] = lambda1[1] + lambda2[1] + lambda1[1]*chi1;
  lsum[2] = lambda1[2] + lambda2[2] + lambda1[2]*chi1;

  double n_eq[19], n_neq[19];

  for (i=0; i<lbmodel.n_veloc; i++) {
    n_eq[i] = w[i] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  				   + scalar(lsum,c[i]) 
  				   + 0.5*SQR(scalar(lambda1,c[i])) );
  }
  
  for (i=0; i<lbmodel.n_veloc; i++) {
    n_neq[i] = lbfluid[0][i][index] - n_eq[i];
  }
 
  //fprintf(stderr,"n_eq=(\n");
  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  fprintf(stderr, "%d %f %f %f\n",i,lbfluid[0][i][index],n_eq[i],n_neq[i]);
  //}
  //fprintf(stderr,")\n");

  int reverse[19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };

  if (lbfields[index].nvec[2] > 0.0) {

    lbfluid[0][5][index]  = n_eq[5]  + n_neq[reverse[5]];
    lbfluid[0][11][index] = n_eq[11] + n_neq[reverse[11]];
    lbfluid[0][14][index] = n_eq[14] + n_neq[reverse[14]];
    lbfluid[0][15][index] = n_eq[15] + n_neq[reverse[15]];
    lbfluid[0][18][index] = n_eq[18] + n_neq[reverse[18]];

  } else {

    lbfluid[0][6][index]  = n_eq[6]  + n_neq[reverse[6]];
    lbfluid[0][12][index] = n_eq[12] + n_neq[reverse[12]];
    lbfluid[0][13][index] = n_eq[13] + n_neq[reverse[13]];
    lbfluid[0][16][index] = n_eq[16] + n_neq[reverse[16]];
    lbfluid[0][17][index] = n_eq[17] + n_neq[reverse[17]];

  }

  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  fprintf(stderr,"%d %f %f\n",i,lbfluid[0][i][index],n_eq[i]);
  //}

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  index_t next[19];
  next[0]  = index;
  next[1]  = index + 1;
  next[2]  = index - 1;
  next[3]  = index + yperiod;
  next[4]  = index - yperiod;
  next[5]  = index + zperiod;
  next[6]  = index - zperiod;
  next[7]  = index + (1 + yperiod);
  next[8]  = index - (1 + yperiod);
  next[9]  = index + (1 - yperiod);
  next[10] = index - (1 - yperiod);
  next[11] = index + (1 + zperiod);
  next[12] = index - (1 + zperiod);
  next[13] = index + (1 - zperiod);
  next[14] = index - (1 - zperiod);
  next[15] = index + (yperiod + zperiod);
  next[16] = index - (yperiod + zperiod);
  next[17] = index + (yperiod - zperiod);
  next[18] = index - (yperiod - zperiod);

  double agrid=lbpar.agrid, tau=lbpar.tau;
  double omega = 2./(6.*lbpar.viscosity*tau/(agrid*agrid)+1.);

  for (i=0; i<lbmodel.n_veloc; i++) {
  
    lbfluid[1][i][next[i]] = lbfluid[0][i][index] - omega * (lbfluid[0][i][index] - n_eq[i]);

  }

  if (lbfields[index].nvec[2] > 0.0) {

    lbfluid[1][6][next[6]] = 0.0;
    lbfluid[1][12][next[12]] = 0.0;
    lbfluid[1][13][next[13]] = 0.0;
    lbfluid[1][16][next[16]] = 0.0;
    lbfluid[1][17][next[17]] = 0.0;

  } else {

    lbfluid[1][5][next[5]] = 0.0;
    lbfluid[1][11][next[11]] = 0.0;
    lbfluid[1][14][next[14]] = 0.0;
    lbfluid[1][15][next[15]] = 0.0;
    lbfluid[1][18][next[18]] = 0.0;

  }

  for (i=0;i<lbmodel.n_veloc;i++) {
    lbfluid[1][i][next[i]] -= lbmodel.coeff[i][0];
  }

  lb_calc_local_j(index, &mode[1]);
//  fprintf(stderr,"j=(%e,%f,%f)\n",mode[1],mode[2],mode[3]);
#endif //if 0
}

MDINLINE void lb_boundary_bb_neq_BGK2(index_t index, double *mode) {
  int i;

  /* \todo Can we do it like this or do we have too use pull scheme? */

  fprintf(stderr,"lb_boundary_bb_neq_BGK2(%+.0f)\n",lbfields[index].nvec[2]);

  /* stream the old ones (they are relaxed already, initial condition?) */
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  index_t next[19];
  next[0]  = index;
  next[1]  = index + 1;
  next[2]  = index - 1;
  next[3]  = index + yperiod;
  next[4]  = index - yperiod;
  next[5]  = index + zperiod;
  next[6]  = index - zperiod;
  next[7]  = index + (1 + yperiod);
  next[8]  = index - (1 + yperiod);
  next[9]  = index + (1 - yperiod);
  next[10] = index - (1 - yperiod);
  next[11] = index + (1 + zperiod);
  next[12] = index - (1 + zperiod);
  next[13] = index + (1 - zperiod);
  next[14] = index - (1 - zperiod);
  next[15] = index + (yperiod + zperiod);
  next[16] = index - (yperiod + zperiod);
  next[17] = index + (yperiod - zperiod);
  next[18] = index - (yperiod - zperiod);

  for (i=0; i<lbmodel.n_veloc; i++) {
    lbfluid[1][i][next[i]] = lbfluid[0][i][index];
  }

#if 0
  das geht so nicht weil beim stroemen nun schon dinge ueberschrieben werden die man auf nachfolgenden Boundary nodes noch braucht!!!
#endif

  /* bb_neq and relax the new ones (leaps a bit ahead then?) */
  mode[0] = lbpar.rho_lb_units;
  for (i=0;i<lbmodel.n_veloc; i++) {
    mode[0] += lbfluid[1][i][index];
  }

  //lb_calc_local_fields(index, &mode[0], &mode[1], NULL);

//  fprintf(stderr, "rho=%f j=(%e,%f,%f)\n",mode[0],mode[1],mode[2],mode[3]);

  mode[1] = 0.0;
  mode[2] = 0.0;
  mode[3] = 0.0;

  for (i=0;i<lbmodel.n_veloc;i++) {
    lbfluid[1][i][index] += lbmodel.coeff[i][0];
  }

  if (lbfields[index].nvec[2] > 0.0) {

    lbfluid[1][5][index] = 0.0;
    lbfluid[1][11][index] = 0.0;
    lbfluid[1][14][index] = 0.0;
    lbfluid[1][15][index] = 0.0;
    lbfluid[1][18][index] = 0.0;

  } else {

    lbfluid[1][6][index] = 0.0;
    lbfluid[1][12][index] = 0.0;
    lbfluid[1][13][index] = 0.0;
    lbfluid[1][16][index] = 0.0;
    lbfluid[1][17][index] = 0.0;

  }

  for (i=0; i<lbmodel.n_veloc; i++) {
    fprintf(stderr,"%d %f\n",i,lbfluid[1][i][index]);
  }

  const double A = lbfields[index].nvec[2] * 1./6.;
  const double P = 5./36.;
  const double Z = 5./108.; // 11./108.; war falsch!

  const double w[19] = { 7./18.,
                         1./12., 1./12., 1./12., 1./12., 1./18., 1./18.,
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36. };

  double (*c)[3] = lbmodel.c;
      
  double chi1, chi2, lambda1[3], lambda2[3], lsum[3];
  double c0;

  lambda1[0] = (mode[1]/mode[0])/lbmodel.c_sound_sq;
  lambda1[1] = (mode[2]/mode[0])/lbmodel.c_sound_sq;
  lambda1[2] = (mode[3]/mode[0])/P;
  
  chi1 = -A*lambda1[2];

  lambda2[0] = 0.0;
  lambda2[1] = 0.0;
  lambda2[2] = -Z/P*SQR(lambda1[2]);

  chi2 = -A*lambda2[2] - 0.5*P*SQR(lambda1[2]) - 1./6.*(SQR(lambda1[0])+SQR(lambda1[1]));

  lsum[0] = lambda1[0] + lambda2[0] + lambda1[0]*chi1;
  lsum[1] = lambda1[1] + lambda2[1] + lambda1[1]*chi1;
  lsum[2] = lambda1[2] + lambda2[2] + lambda1[2]*chi1;

  double n_eq[19], n_neq[19];

  for (i=0; i<lbmodel.n_veloc; i++) {
    n_eq[i] = w[i] * mode[0] * ( 1 + chi1 + chi2 + 0.5*SQR(chi1) 
  				   + scalar(lsum,c[i]) 
  				   + 0.5*SQR(scalar(lambda1,c[i])) );
  }
  
  for (i=0; i<lbmodel.n_veloc; i++) {
    n_neq[i] = lbfluid[0][i][index] - n_eq[i];
  }
 
  //fprintf(stderr,"n_eq=(\n");
  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  fprintf(stderr, "%d %f %f %f\n",i,lbfluid[0][i][index],n_eq[i],n_neq[i]);
  //}
  //fprintf(stderr,")\n");

  int reverse[19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };

  if (lbfields[index].nvec[2] > 0.0) {

    lbfluid[0][5][index]  = n_eq[5]  + n_neq[reverse[5]];
    lbfluid[0][11][index] = n_eq[11] + n_neq[reverse[11]];
    lbfluid[0][14][index] = n_eq[14] + n_neq[reverse[14]];
    lbfluid[0][15][index] = n_eq[15] + n_neq[reverse[15]];
    lbfluid[0][18][index] = n_eq[18] + n_neq[reverse[18]];

  } else {

    lbfluid[0][6][index]  = n_eq[6]  + n_neq[reverse[6]];
    lbfluid[0][12][index] = n_eq[12] + n_neq[reverse[12]];
    lbfluid[0][13][index] = n_eq[13] + n_neq[reverse[13]];
    lbfluid[0][16][index] = n_eq[16] + n_neq[reverse[16]];
    lbfluid[0][17][index] = n_eq[17] + n_neq[reverse[17]];

  }

  for (i=0; i<lbmodel.n_veloc; i++) {
    fprintf(stderr,"%d %f %f\n",i,lbfluid[1][i][index],n_eq[i]);
  }

  double agrid=lbpar.agrid, tau=lbpar.tau;
  double omega = 2./(6.*lbpar.viscosity*tau/(agrid*agrid)+1.);

  for (i=0; i<lbmodel.n_veloc; i++) {
  
    lbfluid[1][i][index] = lbfluid[1][i][index] - omega * (lbfluid[1][i][index] - n_eq[i]);

  }

  if (lbfields[index].nvec[2] > 0.0) {

    lbfluid[1][6][index] = 0.0;
    lbfluid[1][12][index] = 0.0;
    lbfluid[1][13][index] = 0.0;
    lbfluid[1][16][index] = 0.0;
    lbfluid[1][17][index] = 0.0;

  } else {

    lbfluid[1][5][index] = 0.0;
    lbfluid[1][11][index] = 0.0;
    lbfluid[1][14][index] = 0.0;
    lbfluid[1][15][index] = 0.0;
    lbfluid[1][18][index] = 0.0;

  }

  for (i=0;i<lbmodel.n_veloc;i++) {
    lbfluid[1][i][next[i]] -= lbmodel.coeff[i][0];
  }

  //lb_calc_local_j(index, &mode[1]);
  //fprintf(stderr,"j=(%e,%f,%f)\n",mode[1],mode[2],mode[3]);

}

MDINLINE void lb_boundary_relax_modes(index_t index, double *mode, double *pi) {

  double rho, j[3], pi_eq[5], *f = lbfields[index].force;
  double agrid=lbpar.agrid, tau=lbpar.tau;

  double gamma_shear_wall, gamma_bulk_wall, gamma_shear_z;
  double gamma_odd, gamma_even;
  gamma_shear_wall = 1. - 2./(6.*lbpar.viscosity*tau/(agrid*agrid)+1.);
  gamma_bulk_wall = gamma_shear_wall;
  gamma_odd = gamma_even = gamma_shear_wall;
  gamma_shear_z = gamma_shear_wall;

  fprintf(stderr,"lb_boundary_relax_modes()\n");

  rho = mode[0];
  
  /* redefine momentum density for equilibrium parts */
  //j[0] = mode[1] + 0.5*lbfields[index].force[0];
  //j[1] = mode[2] + 0.5*lbfields[index].force[1];
  //j[2] = mode[3] + 0.5*lbfields[index].force[2];

  //if (lbfields[index].nvec[2] < 0.0) {
  //  j[2] = 1./6.*rho;
  //} else {
  //  j[2] = -1./6.*rho;
  //}

  j[0] = mode[1];
  j[1] = mode[2];
  j[2] = mode[3];

  fprintf(stderr,"j=(%f,%f,%f)\n",j[0],j[1],j[2]);

  double A;
  if (lbfields[index].nvec[2] > 0.0) {
    A =  -1./6;
  } else {
    A = 1./6.;
  }

  /* equilibrium part of stress like modes */
  pi_eq[0] = (SQR(j[0])+SQR(j[1]))/rho;
  pi_eq[1] = (SQR(j[0])-SQR(j[1]))/rho;
  pi_eq[2] = j[0]*j[1]/rho;

  pi_eq[3] = j[0]*(j[2]/rho - A);
  pi_eq[4] = j[1]*(j[2]/rho - A);

  fprintf(stderr,"modes=(%f,%f,%f,%f,%e,%e,%e,%e,%e)\n",mode[0],mode[1],mode[2],mode[3],mode[4],mode[5],mode[6],mode[7],mode[8]);
  fprintf(stderr,"pi=(%.6f,%.6f,%.6f,%e,%.6f,%.6f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);
  fprintf(stderr,"pi_eq=(%.9f,%.9f,%.9f,%.9f,%.9f)\n",pi_eq[0],pi_eq[1],pi_eq[2],pi_eq[3],pi_eq[4]);

  /* relax the stress like modes */
  mode[4] = pi_eq[0] + gamma_bulk_wall*(mode[4] - pi_eq[0]);
  mode[5] = pi_eq[1] + gamma_shear_wall*(mode[5] - pi_eq[1]);
  mode[6] = pi_eq[2] + gamma_shear_wall*(mode[6] - pi_eq[2]);

  mode[7] = pi_eq[3] + gamma_shear_z*(mode[7] - pi_eq[3]);
  mode[8] = pi_eq[4] + gamma_shear_z*(mode[8] - pi_eq[4]);

  pi[3] = gamma_shear_wall * pi[3];

  fprintf(stderr,"%e\n",pi[3]/2.);

  //double C = -1./12.*(1.0+gamma_shear_z)*f[0];
  //fprintf(stderr,"C = %e\n",C);
  //mode[7] += C;
  
  /* relax the remaining modes */
  mode[9]  = gamma_odd  * mode[9];
  mode[10] = gamma_odd  * mode[10];
  mode[11] = gamma_odd  * mode[11];
  mode[12] = gamma_odd  * mode[12];
  mode[13] = gamma_even * mode[13];

  fprintf(stderr,"relaxed pi=(%.6f,%.6f,%.6f,%e,%.6f,%.6f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);
  fprintf(stderr,"relaxed modes=(%.9f,%.9f,%.9f,%.9f,%.12f,%.12f,%.9f,%e,%.9f)\n",mode[0],mode[1],mode[2],mode[3],mode[4],mode[5],mode[6],mode[7],mode[8]);

}

MDINLINE void lb_boundary_apply_forces(index_t index, double *mode) {
#if 0 //problem with slip_pref (georg, 03.08.10)
  double *f = lbfields[index].force;
  double rho = mode[0] + lbpar.rho_lb_units;
  double u[3];

  f[0] = lbpar.ext_force[0];
  f[1] = lbpar.ext_force[1];
  f[2] = lbpar.ext_force[2];

  fprintf(stderr,"force = (%e,%f,%f)\n",f[0],f[1],f[2]);
 
  f[0] += -lb_boundary_par.slip_pref * (mode[1]);//+0.5*f[0])/rho;
  f[1] += -lb_boundary_par.slip_pref * (mode[2]);//+0.5*f[1])/rho;
  f[2] += -lb_boundary_par.slip_pref * (mode[3]);//+0.5*f[2])/rho;

  fprintf(stderr,"force = (%e,%f,%f)\n",f[0],f[1],f[2]);
 
  /* redefine velocity field */
  //u[0] = (mode[1] + 0.5*f[0])/rho;
  //u[1] = (mode[2] + 0.5*f[1])/rho;
  //u[2] = (mode[3] + 0.5*f[2])/rho;

  /* correction term */
  //C[0] = (1.+gamma_bulk)*u[0]*f[0] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  //C[2] = (1.+gamma_bulk)*u[1]*f[1] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  //C[5] = (1.+gamma_bulk)*u[2]*f[2] + 1./3.*(gamma_bulk-gamma_shear)*scalar(u,f);
  //C[1] = 1./2.*(1.+gamma_shear)*(u[0]*f[1]+u[1]*f[0]);
  //C[3] = 1./2.*(1.+gamma_shear)*(u[0]*f[2]+u[2]*f[0]);
  //C[4] = 1./2.*(1.+gamma_shear)*(u[1]*f[2]+u[2]*f[1]);

  /* update momentum modes */
  mode[1] += f[0];
  mode[2] += f[1];
  mode[3] += f[2];

  /* update stress modes */
  //mode[4] += C[0] + C[2] + C[5];
  //mode[5] += 2.*C[0] - C[2] - C[5];
  //mode[6] += C[2] - C[5];
  //mode[7] += C[1];
  //mode[8] += C[3];
  //mode[9] += C[4];
#endif //if 0
}

MDINLINE void lb_boundary_calc_n_push(index_t index, double *m) {
  int i;

  int n_wall_veloc = 14;
  double w[19]    = { 7./18., 
                      1./12., 1./12., 1./12., 1./12., 1./18., 1./18.,
                      1./36., 1./36., 1./36., 1./36., 1./36., 1./36.,
                      1./36., 1./36., 1./36., 1./36., 1./36., 1./36. };
  double norm[14] = { 1., 1./3., 1./3., 5./36., 
                      4./9., 4./9., 1./9., 5./108., 5./108., 
		      1./15., 1./15., 11./324., 1./12., 28./165. };

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  index_t next[19];
  next[0]  = index;
  next[1]  = index + 1;
  next[2]  = index - 1;
  next[3]  = index + yperiod;
  next[4]  = index - yperiod;
  next[5]  = index + zperiod;
  next[6]  = index - zperiod;
  next[7]  = index + (1 + yperiod);
  next[8]  = index - (1 + yperiod);
  next[9]  = index + (1 - yperiod);
  next[10] = index - (1 - yperiod);
  next[11] = index + (1 + zperiod);
  next[12] = index - (1 + zperiod);
  next[13] = index + (1 - zperiod);
  next[14] = index - (1 - zperiod);
  next[15] = index + (yperiod + zperiod);
  next[16] = index - (yperiod + zperiod);
  next[17] = index + (yperiod - zperiod);
  next[18] = index - (yperiod - zperiod);

  fprintf(stderr,"lb_boundary_calc_n_push()\n");
  
  //fprintf(stderr,"(%f,%e,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) %f\n",m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12],m[13],lbfields[index].nvec[2]);


  for (i=0;i<n_wall_veloc;i++) {
    m[i] = m[i]/norm[i];
  }

  m[3] *= 1./6.;
  m[4] *= 1./3.;
  m[7] *= 1./6.;
  m[8] *= 1./6.;
  m[9] *= 1./5.;
  m[10] *= 1./5.;
  m[11] *= 1./36.;
  m[12] *= 1./4.;
  m[13] *= 1./55.;

  if (lbfields[index].nvec[2] > 0.0) {

    lbfluid[1][0][next[0]] = m[0] + m[3] - 2.*m[4] - 2.*m[11] - 12.*m[13];
    lbfluid[1][1][next[1]] = m[0] + m[1] + m[3] + m[4] + m[5] + m[7] - 2.*m[9] + m[11] + m[12] - 28.*m[13];
    lbfluid[1][2][next[2]] = m[0] - m[1] + m[3] + m[4] + m[5] - m[7] + 2.*m[9] + m[11] + m[12] - 28.*m[13];
    lbfluid[1][3][next[3]] = m[0] + m[2] + m[3] + m[4] - m[5] + m[8] - 2.*m[10] + m[11] - m[12] - 28.*m[13];
    lbfluid[1][4][next[4]] = m[0] - m[2] + m[3] + m[4] - m[5] - m[8] + 2.*m[10] + m[11] - m[12] - 28.*m[13];

    lbfluid[1][7][next[7]] = m[0] + m[1] + m[2] + m[3] + 4.*m[4] + m[6] + m[7] + m[8] + 3.*m[9] + 3.*m[10] + 4.*m[11] + 42.*m[13];
    lbfluid[1][8][next[8]] = m[0] - m[1] - m[2] + m[3] + 4.*m[4] + m[6] - m[7] - m[8] - 3.*m[9] - 3.*m[10] + 4.*m[11] + 42.*m[13];
    lbfluid[1][9][next[9]] = m[0] + m[1] - m[2] + m[3] + 4.*m[4] - m[6] + m[7] - m[8] + 3.*m[9] - 3.*m[10] + 4.*m[11] + 42.*m[13];
    lbfluid[1][10][next[10]] = m[0] - m[1] + m[2] + m[3] + 4.*m[4] - m[6] - m[7] + m[8] - 3.*m[9] + 3.*m[10] + 4.*m[11] + 42.*m[13];

    lbfluid[1][5][next[5]] = m[0] - 5.*m[3] - 2.*m[4] + 22.*m[11]; 

    lbfluid[1][11][next[11]] = m[0] + m[1] - 5.*m[3] + m[4] + m[5] - 5.*m[7] - 11.*m[11] - 3.*m[12];
    lbfluid[1][14][next[14]] = m[0] - m[1] - 5.*m[3] + m[4] + m[5] + 5.*m[7] - 11.*m[11] - 3.*m[12];
    lbfluid[1][15][next[15]] = m[0] + m[2] - 5.*m[3] + m[4] - m[5] - 5.*m[8] - 11.*m[11] + 3.*m[12];
    lbfluid[1][18][next[18]] = m[0] - m[2] - 5.*m[3] + m[4] - m[5] + 5.*m[8] - 11.*m[11] + 3.*m[12];
 
    lbfluid[1][6][next[6]] = 0.0;
    lbfluid[1][12][next[12]] = 0.0;
    lbfluid[1][13][next[13]] = 0.0;
    lbfluid[1][16][next[16]] = 0.0;
    lbfluid[1][17][next[17]] = 0.0;

  } else {

    lbfluid[1][0][next[0]] = m[0] - m[3] - 2.*m[4] + 2.*m[11] + 12.*m[13];
    lbfluid[1][1][next[1]] = m[0] + m[1] - m[3] + m[4] + m[5] - m[7] - 2.*m[9] - m[11] - m[12] - 28.*m[13];
    lbfluid[1][2][next[2]] = m[0] - m[1] - m[3] + m[4] + m[5] + m[7] + 2.*m[9] - m[11] - m[12] - 28.*m[13];
    lbfluid[1][3][next[3]] = m[0] + m[2] - m[3] + m[4] - m[5] - m[8] - 2.*m[10] - m[11] + m[12] - 28.*m[13];
    lbfluid[1][4][next[4]] = m[0] - m[2] - m[3] + m[4] - m[5] + m[8] + 2.*m[10] - m[11] + m[12] - 28.*m[13];

    lbfluid[1][7][next[7]] = m[0] + m[1] + m[2] - m[3] + 4.*m[4] + m[6] - m[7] - m[8] + 3.*m[9] + 3.*m[10] - 4.*m[11] + 42.*m[13];
    lbfluid[1][8][next[8]] = m[0] - m[1] - m[2] - m[3] + 4.*m[4] + m[6] + m[7] + m[8] - 3.*m[9] - 3.*m[10] - 4.*m[11] + 42.*m[13];
    lbfluid[1][9][next[9]] = m[0] + m[1] - m[2] - m[3] + 4.*m[4] - m[6] - m[7] + m[8] + 3.*m[9] - 3.*m[10] - 4.*m[11] + 42.*m[13];
    lbfluid[1][10][next[10]] = m[0] - m[1] + m[2] - m[3] + 4.*m[4] - m[6] + m[7] - m[8] - 3.*m[9] + 3.*m[10] - 4.*m[11] + 42.*m[13];

    lbfluid[1][6][next[6]] = m[0] + 5.*m[3] - 2.*m[4] - 22.*m[11]; 

    lbfluid[1][13][next[13]] = m[0] + m[1] + 5.*m[3] + m[4] + m[5] + 5.*m[7] + 11.*m[11] + 3.*m[12];
    lbfluid[1][12][next[12]] = m[0] - m[1] + 5.*m[3] + m[4] + m[5] - 5.*m[7] + 11.*m[11] + 3.*m[12];
    lbfluid[1][17][next[17]] = m[0] + m[2] + 5.*m[3] + m[4] - m[5] + 5.*m[8] + 11.*m[11] - 3.*m[12];
    lbfluid[1][16][next[16]] = m[0] - m[2] + 5.*m[3] + m[4] - m[5] - 5.*m[8] + 11.*m[11] - 3.*m[12];
 
    lbfluid[1][5][next[5]] = 0.0;
    lbfluid[1][11][next[11]] = 0.0;
    lbfluid[1][14][next[14]] = 0.0;
    lbfluid[1][15][next[15]] = 0.0;
    lbfluid[1][18][next[18]] = 0.0;

  }

  //fprintf(stderr,"( ");
  for (i=0;i<lbmodel.n_veloc;i++) {
    lbfluid[1][i][next[i]] *= w[i];
  //  fprintf(stderr,"%f ",lbfluid[1][i][next[i]]);
    lbfluid[1][i][next[i]] -= lbmodel.coeff[i][0];
  }
  //fprintf(stderr,") %f\n",lbfields[index].nvec[2]);

  //double **tmp = lbfluid[0];
  //lbfluid[0] = lbfluid[1];
  //double rho, j[3], agrid=lbpar.agrid, tau=lbpar.tau, *f=lbfields[index].force;
  //
  //lb_calc_local_fields(index, &rho, j, NULL, 0);
  ////lb_calc_local_rho(index, &rho);
  ////lb_calc_local_j(index,j);
  //lbfluid[0] = tmp;
  //fprintf(stderr,"%p %d (%.12e,%f,%f) %.12e\n",lbfluid[1],index,(j[0]+f[0])/rho,(j[1]+f[1])/rho,(j[2]+f[2])/rho,rho);

}

MDINLINE void lb_boundary_equilibrium_push(int index, double *mode) {
  int i;

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[19];
  next[0]  = index;
  next[1]  = index + 1;
  next[2]  = index - 1;
  next[3]  = index + yperiod;
  next[4]  = index - yperiod;
  next[5]  = index + zperiod;
  next[6]  = index - zperiod;
  next[7]  = index + (1 + yperiod);
  next[8]  = index - (1 + yperiod);
  next[9]  = index + (1 - yperiod);
  next[10] = index - (1 - yperiod);
  next[11] = index + (1 + zperiod);
  next[12] = index - (1 + zperiod);
  next[13] = index + (1 - zperiod);
  next[14] = index - (1 - zperiod);
  next[15] = index + (yperiod + zperiod);
  next[16] = index - (yperiod + zperiod);
  next[17] = index + (yperiod - zperiod);
  next[18] = index - (yperiod - zperiod);
  
#ifdef D3Q19
  const double A = 1./6.;
  const double P = 5./36.;
  const double Z = 5./108.;  // 11./108.; war falsch!

  const double w[19] = { 7./18.,
                         1./12., 1./12., 1./12., 1./12., 1./18., 1./18.,
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36., 
			 1./36., 1./36., 1./36., 1./36. };

  double (*c)[3] = lbmodel.c;

  double *nvec = lbfields[index].nvec;
      
  double chi1, chi2, lambda1[3], lambda2[3], lsum[3];
  double c0;

  double rho = mode[0] + lbpar.rho_lb_units;

  lambda1[0] = (mode[1]+A*nvec[0]) / (lbmodel.c_sound_sq*rho);
  lambda1[1] = (mode[2]+A*nvec[1]) / (lbmodel.c_sound_sq*rho);
  lambda1[2] = (mode[3]+A*nvec[2]) / (lbmodel.c_sound_sq*rho);

  chi1 = -A*lambda1[2];

  lambda2[0] = 0.0;
  lambda2[1] = 0.0;
  lambda2[2] = -Z/P*SQR(lambda1[2]);

  chi2 = -A*lambda2[2] - 0.5*P*SQR(lambda1[2]) - 1./6.*(SQR(lambda1[0])+SQR(lambda1[1]));

  c0 = chi1 + 0.5*SQR(chi1) + chi2;

  lsum[0] = lambda1[0] + lambda2[0] + lambda1[0]*chi1;
  lsum[1] = lambda1[1] + lambda2[1] + lambda1[1]*chi1;
  lsum[2] = lambda1[2] + lambda2[2] + lambda1[2]*chi1;

  //if ((index == get_linear_index(1,1,21,lblattice.halo_grid))
  //     || (index == get_linear_index(1,1,1,lblattice.halo_grid))) {
  //  fprintf(stderr, "rho=%e j=(%e,%e,%e)\n",rho,mode[1],mode[2],mode[3]);
  //  fprintf(stderr,"chi1 = %e chi2 = %f lambda1 = (%e,%e,%e) lambda2 = (%e,%e,%e)\n",chi1,chi2,lambda1[0],lambda1[1],lambda1[2],lambda2[0],lambda2[1],lambda2[2]);
  //}

  for (i=0; i<lbmodel.n_veloc; i++) {

    if (scalar(c[i],lbfields[index].nvec) >= 0.0) {
      
      lbfluid[1][i][next[i]] = w[i] * rho
			     + w[i] * rho * ( c0 
					      + scalar(lsum,c[i]) 
					      + 0.5*SQR(scalar(lambda1,c[i])) )
	- lbmodel.coeff[i][0]*lbpar.rho_lb_units;

      if ((index == get_linear_index(1,1,21,lblattice.halo_grid))
	   || (index == get_linear_index(1,1,1,lblattice.halo_grid))) {
	//fprintf(stderr,"[%d] %e\n",i,lbfluid[1][i][next[i]]+lbmodel.coeff[i][0]*lbpar.rho_lb_units);
      }

    } else {

      /* nothing streams through the wall 
       * this is important for the halo regions */
      lbfluid[1][i][next[i]] = -lbmodel.coeff[i][0]*lbpar.rho_lb_units;

    }

  }
#endif
}

MDINLINE void lb_boundary_collisions(int index, double *modes) {

  double rho, v[3], f[3];

  double pi[6];

  lb_boundary_bb_neq_BGK(index, modes);
  return;

  //lb_boundary_calc_modes(index, modes, pi);

  //rho = modes[0];// + lbpar.rho_lb_units;
  //v[0] = modes[1]/rho;
  //v[1] = modes[2]/rho;
  //v[2] = modes[3]/rho;
  //
  //fprintf(stderr,"%f\n",lb_boundary_par.slip_pref);
  //
  //f[0] = - lb_boundary_par.slip_pref*v[0];
  //f[1] = 0.0;
  //f[2] = 0.0;

  //lbfields[index].force[0] += f[0];
  //lbfields[index].force[1] += f[1];
  //lbfields[index].force[2] += f[2];
  //lbfields[index].force[0] = -modes[1];
  //lbfields[index].force[1] = -modes[2];
  //lbfields[index].force[2] = 0.0;

  //modes[1] = 0.0;
  //modes[2] = 0.0;

  //lb_boundary_relax_modes(index, modes, pi);

  //lb_boundary_apply_forces(index, modes);

  //modes[1] = 0.0;
  //modes[2] = 0.0;

  //lb_boundary_calc_n_push(index, modes);

}

MDINLINE void lb_partial_slip() {
  index_t index;
  double modes[19];

  for (index=0; index<lblattice.halo_grid_volume; index++) {

    /* improve on this loop with lists - How does Tony do it? */

    if (lbfields[index].boundary) {

      lb_boundary_bb_neq_BGK(index,modes);

    }

  }

}

#if 0
#if 0
/* leaves populations in the wall */
MDINLINE void lb_bounce_back2() {

#ifdef D3Q19
  int k;
  int reverse[] = { 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16 };
  double new_n[18];

  for (k=0;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      /* bounce back to lower indices */
      new_n[reverse[0]]  = lbfluid[k].n[0];
      new_n[reverse[2]]  = lbfluid[k].n[2];
      new_n[reverse[4]]  = lbfluid[k].n[4];
      new_n[reverse[6]]  = lbfluid[k].n[6];
      new_n[reverse[9]]  = lbfluid[k].n[9];
      new_n[reverse[10]] = lbfluid[k].n[10];
      new_n[reverse[13]] = lbfluid[k].n[13];
      new_n[reverse[14]] = lbfluid[k].n[14];
      new_n[reverse[17]] = lbfluid[k].n[17];
      
      /* bounce back to higher indices */
      new_n[reverse[1]]  = lbfluid[k].n[1];
      new_n[reverse[3]]  = lbfluid[k].n[3];
      new_n[reverse[5]]  = lbfluid[k].n[5];
      new_n[reverse[7]]  = lbfluid[k].n[7];
      new_n[reverse[8]]  = lbfluid[k].n[8];
      new_n[reverse[11]] = lbfluid[k].n[11];
      new_n[reverse[12]] = lbfluid[k].n[12];
      new_n[reverse[15]] = lbfluid[k].n[15];
      new_n[reverse[16]] = lbfluid[k].n[16];

      memcpy(lbfluid[k].n,new_n,18*sizeof(double));
      
    }

  }
#else
#error Bounce back boundary conditions are only implemented for D3Q18!
#endif
}
#endif

#define CANONICAL
#define NONE
MDINLINE void lb_partial_slip() {

  int i,k;

  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   0;                       // ( 1, 0, 0) +
  next[1]  = - 0;                       // (-1, 0, 0)
  next[2]  =   0;                       // ( 0, 1, 0) +
  next[3]  = - 0;                       // ( 0,-1, 0)
  next[4]  =   zperiod;                 // ( 0, 0, 1) +
  next[5]  = - zperiod;                 // ( 0, 0,-1)
  next[6]  =   0;                       // ( 1, 1, 0) +
  next[7]  = - 0;                       // (-1,-1, 0)
  next[8]  =   0;                       // ( 1,-1, 0) 
  next[9]  = - 0;                       // (-1, 1, 0) +
  next[10] =   zperiod;                 // ( 1, 0, 1) +
  next[11] = - zperiod;                 // (-1, 0,-1)
  next[12] = - zperiod;                 // ( 1, 0,-1)
  next[13] =   zperiod;                 // (-1, 0, 1) +
  next[14] =   zperiod;                 // ( 0, 1, 1) +
  next[15] = - zperiod;                 // ( 0,-1,-1)
  next[16] = - zperiod;                 // ( 0, 1,-1)
  next[17] =   zperiod;                 // ( 0,-1, 1) +

  int reflect[] = { 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 16, 17, 14, 15 };

  double delta_j[3], u[3], u1[3], u2[3], u3[3];
  double gamma = lb_boundary_par.slip_pref;

  for (k=0;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary == -1) {

      /* calculate fields at nodes used for extrapolation */
      lb_calc_local_fields(&lbfluid[k-zperiod],0);
      lb_calc_local_fields(&lbfluid[k-2*zperiod],0);
      lb_calc_local_fields(&lbfluid[k-3*zperiod],0);

      u1[0] = lbfluid[k-3*zperiod].j[0] / *lbfluid[k-3*zperiod].rho;
      u1[1] = lbfluid[k-3*zperiod].j[1] / *lbfluid[k-3*zperiod].rho;
      u2[0] = lbfluid[k-2*zperiod].j[0] / *lbfluid[k-2*zperiod].rho;
      u2[1] = lbfluid[k-2*zperiod].j[1] / *lbfluid[k-2*zperiod].rho;
      u3[0] = lbfluid[k-zperiod].j[0] / *lbfluid[k-zperiod].rho;
      u3[1] = lbfluid[k-zperiod].j[1] / *lbfluid[k-zperiod].rho;

#ifdef NONE
      u[0] = u3[0];
      u[1] = u3[1];
#endif
#ifdef LINEAR
      u[0] = u3[0] + 0.5*(u3[0]-u2[0]);
      u[1] = u3[1] + 0.5*(u3[1]-u2[1]);
#endif
#ifdef QUADRATIC
      /* extrapolation for velocity at wall */
      u[0] = 15./8.*u3[0] - 5./4.*u2[0] + 3./8.*u1[0];
      u[1] = 15./8.*u3[1] - 5./4.*u2[1] + 3./8.*u1[1];
#endif

      /* calculate momentum transfer due to friction */
      delta_j[0] = - gamma * u[0] * lbpar.tau;
      delta_j[1] = - gamma * u[1] * lbpar.tau;
  
      /* apply friction force */
#ifdef CANONICAL
      lbfluid[k-zperiod].n[0]  +=   1./6.*delta_j[0];
      lbfluid[k-zperiod].n[1]  += - 1./6.*delta_j[0];
      lbfluid[k-zperiod].n[2]  +=   1./6.*delta_j[1];
      lbfluid[k-zperiod].n[3]  += - 1./6.*delta_j[1];
      lbfluid[k-zperiod].n[4]  +=   0.0;
      lbfluid[k-zperiod].n[5]  +=   0.0;
      lbfluid[k-zperiod].n[6]  +=   1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k-zperiod].n[7]  += - 1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k-zperiod].n[8]  +=   1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k-zperiod].n[9]  += - 1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k-zperiod].n[10] +=   1./12.*delta_j[0];
      lbfluid[k-zperiod].n[11] += - 1./12.*delta_j[0];
      lbfluid[k-zperiod].n[12] +=   1./12.*delta_j[0];
      lbfluid[k-zperiod].n[13] += - 1./12.*delta_j[0];
      lbfluid[k-zperiod].n[14] +=   1./12.*delta_j[1];
      lbfluid[k-zperiod].n[15] += - 1./12.*delta_j[1];
      lbfluid[k-zperiod].n[16] +=   1./12.*delta_j[1];
      lbfluid[k-zperiod].n[17] += - 1./12.*delta_j[1];          
#else
      lbfluid[k-zperiod].n[0]  +=   0.0;//1./6.*delta_j[0];
      lbfluid[k-zperiod].n[1]  += - 0.0;//1./6.*delta_j[0];
      lbfluid[k-zperiod].n[2]  +=   0.0;//1./6.*delta_j[1];
      lbfluid[k-zperiod].n[3]  += - 0.0;//1./6.*delta_j[1];
      lbfluid[k-zperiod].n[4]  +=   0.0;
      lbfluid[k-zperiod].n[5]  +=   0.0;
      lbfluid[k-zperiod].n[6]  +=   0.0;//1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k-zperiod].n[7]  += - 0.0;//1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k-zperiod].n[8]  +=   0.0;//1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k-zperiod].n[9]  += - 0.0;//1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k-zperiod].n[10] +=   1./2.*delta_j[0];
      lbfluid[k-zperiod].n[11] += - 0.0;//1./12.*delta_j[0];
      lbfluid[k-zperiod].n[12] +=   0.0;//1./12.*delta_j[0];
      lbfluid[k-zperiod].n[13] += - 1./2.*delta_j[0];
      lbfluid[k-zperiod].n[14] +=   1./2.*delta_j[1];
      lbfluid[k-zperiod].n[15] += - 0.0;//1./12.*delta_j[1];
      lbfluid[k-zperiod].n[16] +=   0.0;//1./12.*delta_j[1];
      lbfluid[k-zperiod].n[17] += - 1./2.*delta_j[1];          
#endif
      
      for (i=0;i<18;i++) {
	if (lbfluid[k-zperiod].n[i] < 0.0) {
	  fprintf(stderr,"negative %d\n",i);
	}
      }
      
     }

    if (lbfluid[k].boundary == 1) {

      /* calculate fields at nodes used for extrapolation */
      lb_calc_local_fields(&lbfluid[k+zperiod],0);
      lb_calc_local_fields(&lbfluid[k+2*zperiod],0);
      lb_calc_local_fields(&lbfluid[k+3*zperiod],0);

      u1[0] = lbfluid[k+3*zperiod].j[0] / *lbfluid[k+3*zperiod].rho;
      u1[1] = lbfluid[k+3*zperiod].j[1] / *lbfluid[k+3*zperiod].rho;
      u2[0] = lbfluid[k+2*zperiod].j[0] / *lbfluid[k+2*zperiod].rho;
      u2[1] = lbfluid[k+2*zperiod].j[1] / *lbfluid[k+2*zperiod].rho;
      u3[0] = lbfluid[k+zperiod].j[0] / *lbfluid[k+zperiod].rho;
      u3[1] = lbfluid[k+zperiod].j[1] / *lbfluid[k+zperiod].rho;

#ifdef NONE
      u[0] = u3[0];
      u[1] = u3[1];
#endif
#ifdef LINEAR
      u[0] = u3[0] - 0.5*(u2[0]-u3[0]);
      u[1] = u3[1] - 0.5*(u2[1]-u3[1]);
#endif
#ifdef QUADRATIC
      /* extrapolation for velocity at wall */
      u[0] = 15./8.*u3[0] - 5./4.*u2[0] + 3./8.*u1[0];
      u[1] = 15./8.*u3[1] - 5./4.*u2[1] + 3./8.*u1[1];
#endif

      /* calculate momentum transfer due to friction */
      delta_j[0] = - gamma * u[0] * lbpar.tau;
      delta_j[1] = - gamma * u[1] * lbpar.tau;
  
      /* apply friction force */
#ifdef CANONICAL
      lbfluid[k+zperiod].n[0]  +=   1./6.*delta_j[0];
      lbfluid[k+zperiod].n[1]  += - 1./6.*delta_j[0];
      lbfluid[k+zperiod].n[2]  +=   1./6.*delta_j[1];
      lbfluid[k+zperiod].n[3]  += - 1./6.*delta_j[1];
      lbfluid[k+zperiod].n[4]  +=   0.0;
      lbfluid[k+zperiod].n[5]  +=   0.0;
      lbfluid[k+zperiod].n[6]  +=   1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k+zperiod].n[7]  += - 1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k+zperiod].n[8]  +=   1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k+zperiod].n[9]  += - 1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k+zperiod].n[10] +=   1./12.*delta_j[0];
      lbfluid[k+zperiod].n[11] += - 1./12.*delta_j[0];
      lbfluid[k+zperiod].n[12] +=   1./12.*delta_j[0];
      lbfluid[k+zperiod].n[13] += - 1./12.*delta_j[0];
      lbfluid[k+zperiod].n[14] +=   1./12.*delta_j[1];
      lbfluid[k+zperiod].n[15] += - 1./12.*delta_j[1];
      lbfluid[k+zperiod].n[16] +=   1./12.*delta_j[1];
      lbfluid[k+zperiod].n[17] += - 1./12.*delta_j[1];          
#else
      lbfluid[k+zperiod].n[0]  +=   0.0;//1./6.*delta_j[0];
      lbfluid[k+zperiod].n[1]  += - 0.0;//1./6.*delta_j[0];
      lbfluid[k+zperiod].n[2]  +=   0.0;//1./6.*delta_j[1];
      lbfluid[k+zperiod].n[3]  += - 0.0;//1./6.*delta_j[1];
      lbfluid[k+zperiod].n[4]  +=   0.0;
      lbfluid[k+zperiod].n[5]  +=   0.0;
      lbfluid[k+zperiod].n[6]  +=   0.0;//1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k+zperiod].n[7]  += - 0.0;//1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k+zperiod].n[8]  +=   0.0;//1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k+zperiod].n[9]  += - 0.0;//1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k+zperiod].n[10] +=   0.0;//1./12.*delta_j[0];
      lbfluid[k+zperiod].n[11] += - 1./2.*delta_j[0];
      lbfluid[k+zperiod].n[12] +=   1./2.*delta_j[0];
      lbfluid[k+zperiod].n[13] += - 0.0;//1./12.*delta_j[0];
      lbfluid[k+zperiod].n[14] +=   0.0;//1./12.*delta_j[1];
      lbfluid[k+zperiod].n[15] += - 1./2.*delta_j[1];
      lbfluid[k+zperiod].n[16] +=   1./2.*delta_j[1];
      lbfluid[k+zperiod].n[17] += - 0.0;//1./12.*delta_j[1];          
#endif

      for (i=0;i<18;i++) {
	if (lbfluid[k+zperiod].n[i] < 0.0) {
	  fprintf(stderr,"negative %d\n",i);
	}
      }
      
    }

  }

  for (k=zperiod;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary == -1) {

      /* reflect from lower indices */
      lbfluid[k].n[reflect[4]]  = lbfluid[k-next[4]].n[4];
      lbfluid[k].n[reflect[10]] = lbfluid[k-next[10]].n[10];
      lbfluid[k].n[reflect[13]] = lbfluid[k-next[13]].n[13];
      lbfluid[k].n[reflect[14]] = lbfluid[k-next[14]].n[14];
      lbfluid[k].n[reflect[17]] = lbfluid[k-next[17]].n[17];

      //lbfluid[k].n[reflect[4]]  = 0.0;
      //lbfluid[k].n[reflect[10]] = 0.0;
      //lbfluid[k].n[reflect[13]] = 0.0;
      //lbfluid[k].n[reflect[14]] = 0.0;
      //lbfluid[k].n[reflect[17]] = 0.0;

      lbfluid[k-next[4]].n[4]   = 0.0;
      lbfluid[k-next[10]].n[10] = 0.0;
      lbfluid[k-next[13]].n[13] = 0.0;
      lbfluid[k-next[14]].n[14] = 0.0;
      lbfluid[k-next[17]].n[17] = 0.0;
      
    }
    
  }

  for (k=lblattice.halo_grid_volume-1-zperiod;k>=0;k--) {

    if (lbfluid[k].boundary == 1) {

      /* reflect from higher indices */
      lbfluid[k].n[reflect[5]]  = lbfluid[k-next[5]].n[5];
      lbfluid[k].n[reflect[11]] = lbfluid[k-next[11]].n[11];
      lbfluid[k].n[reflect[12]] = lbfluid[k-next[12]].n[12];
      lbfluid[k].n[reflect[15]] = lbfluid[k-next[15]].n[15];
      lbfluid[k].n[reflect[16]] = lbfluid[k-next[16]].n[16];

      //lbfluid[k].n[reflect[5]]  = 0.0;
      //lbfluid[k].n[reflect[11]] = 0.0;
      //lbfluid[k].n[reflect[12]] = 0.0;
      //lbfluid[k].n[reflect[15]] = 0.0;
      //lbfluid[k].n[reflect[16]] = 0.0;

      lbfluid[k-next[5]].n[5]   = 0.0;
      lbfluid[k-next[11]].n[11] = 0.0;
      lbfluid[k-next[12]].n[12] = 0.0;
      lbfluid[k-next[15]].n[15] = 0.0;
      lbfluid[k-next[16]].n[16] = 0.0;      
    
    }

  }

}




#if 0
/* friction is applied prior to bouncing */
MDINLINE void lb_slip_boundaries6() {

  int k;
  double *local_n, local_rho, local_j[3];
  double delta_j[3];
  double gamma = lb_boundary_par.slip_pref;

  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   0;                       // ( 1, 0, 0) +
  next[1]  = - 0;                       // (-1, 0, 0)
  next[2]  =   0;                       // ( 0, 1, 0) +
  next[3]  = - 0;                       // ( 0,-1, 0)
  next[4]  =   zperiod;                 // ( 0, 0, 1) +
  next[5]  = - zperiod;                 // ( 0, 0,-1)
  next[6]  =   0;                       // ( 1, 1, 0) +
  next[7]  = - 0;                       // (-1,-1, 0)
  next[8]  =   0;                       // ( 1,-1, 0) 
  next[9]  = - 0;                       // (-1, 1, 0) +
  next[10] =   zperiod;                 // ( 1, 0, 1) +
  next[11] = - zperiod;                 // (-1, 0,-1)
  next[12] = - zperiod;                 // ( 1, 0,-1)
  next[13] =   zperiod;                 // (-1, 0, 1) +
  next[14] =   zperiod;                 // ( 0, 1, 1) +
  next[15] = - zperiod;                 // ( 0,-1,-1)
  next[16] = - zperiod;                 // ( 0, 1,-1)
  next[17] =   zperiod;                 // ( 0,-1, 1) +

  int reflect[] = { 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 16, 17, 14, 15 };

  for (k=0;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary == -1) {

      local_n = lbfluid[k-zperiod].n;

      local_rho = local_n[4] + local_n[10] + local_n[13] + local_n[14] + local_n[17];

      if (local_rho > ROUND_ERROR_PREC) {

	local_j[0] = local_n[10] - local_n[13];
	local_j[1] = local_n[14] - local_n[17];

	/* calculate friction at the wall (lattice units) */
	delta_j[0] = - gamma * local_j[0]/(local_rho) * lbpar.tau;
	delta_j[1] = - gamma * local_j[1]/(local_rho) * lbpar.tau;

	/* apply friction force */
	lbfluid[k-zperiod].n[0]  +=   1./6.*delta_j[0];
	lbfluid[k-zperiod].n[1]  += - 1./6.*delta_j[0];
	lbfluid[k-zperiod].n[2]  +=   1./6.*delta_j[1];
	lbfluid[k-zperiod].n[3]  += - 1./6.*delta_j[1];
	lbfluid[k-zperiod].n[4]  +=   0.0;
	lbfluid[k-zperiod].n[5]  +=   0.0;
	lbfluid[k-zperiod].n[6]  +=   1./12.*(delta_j[0]+delta_j[1]);
	lbfluid[k-zperiod].n[7]  += - 1./12.*(delta_j[0]+delta_j[1]);
	lbfluid[k-zperiod].n[8]  +=   1./12.*(delta_j[0]-delta_j[1]);
	lbfluid[k-zperiod].n[9]  += - 1./12.*(delta_j[0]-delta_j[1]);
	lbfluid[k-zperiod].n[10] +=   1./12.*delta_j[0];
	lbfluid[k-zperiod].n[11] += - 1./12.*delta_j[0];
	lbfluid[k-zperiod].n[12] +=   1./12.*delta_j[0];
	lbfluid[k-zperiod].n[13] += - 1./12.*delta_j[0];
	lbfluid[k-zperiod].n[14] +=   1./12.*delta_j[1];
	lbfluid[k-zperiod].n[15] += - 1./12.*delta_j[1];
	lbfluid[k-zperiod].n[16] +=   1./12.*delta_j[1];
	lbfluid[k-zperiod].n[17] += - 1./12.*delta_j[1];     

      }

    }

    else if (lbfluid[k].boundary == 1) {
      
      local_n = lbfluid[k+zperiod].n;

      local_rho = local_n[5] + local_n[11] + local_n[12] + local_n[15] + local_n[16];

      if (local_rho > ROUND_ERROR_PREC) {

	local_j[0] = - local_n[11] + local_n[12];
	local_j[1] = - local_n[15] + local_n[16];

	/* calculate friction at the wall (lattice units) */
	delta_j[0] = - gamma * local_j[0]/(local_rho) * lbpar.tau;
	delta_j[1] = - gamma * local_j[1]/(local_rho) * lbpar.tau;

	/* apply friction force */
	lbfluid[k+zperiod].n[0]  +=   1./6.*delta_j[0];
	lbfluid[k+zperiod].n[1]  += - 1./6.*delta_j[0];
	lbfluid[k+zperiod].n[2]  +=   1./6.*delta_j[1];
	lbfluid[k+zperiod].n[3]  += - 1./6.*delta_j[1];
	lbfluid[k+zperiod].n[4]  +=   0.0;
	lbfluid[k+zperiod].n[5]  +=   0.0;
	lbfluid[k+zperiod].n[6]  +=   1./12.*(delta_j[0]+delta_j[1]);
	lbfluid[k+zperiod].n[7]  += - 1./12.*(delta_j[0]+delta_j[1]);
	lbfluid[k+zperiod].n[8]  +=   1./12.*(delta_j[0]-delta_j[1]);
	lbfluid[k+zperiod].n[9]  += - 1./12.*(delta_j[0]-delta_j[1]);
	lbfluid[k+zperiod].n[10] +=   1./12.*delta_j[0];
	lbfluid[k+zperiod].n[11] += - 1./12.*delta_j[0];
	lbfluid[k+zperiod].n[12] +=   1./12.*delta_j[0];
	lbfluid[k+zperiod].n[13] += - 1./12.*delta_j[0];
	lbfluid[k+zperiod].n[14] +=   1./12.*delta_j[1];
	lbfluid[k+zperiod].n[15] += - 1./12.*delta_j[1];
	lbfluid[k+zperiod].n[16] +=   1./12.*delta_j[1];
	lbfluid[k+zperiod].n[17] += - 1./12.*delta_j[1];     

      }

    }

  }

  /* bottom-up sweep */
  for (k=zperiod;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary == -1) {

      /* reflect upgoing velocities from lower indices */
      lbfluid[k].n[reflect[4]]  = lbfluid[k-next[4]].n[4];
      lbfluid[k].n[reflect[10]] = lbfluid[k-next[10]].n[10];
      lbfluid[k].n[reflect[13]] = lbfluid[k-next[13]].n[13];
      lbfluid[k].n[reflect[14]] = lbfluid[k-next[14]].n[14];
      lbfluid[k].n[reflect[17]] = lbfluid[k-next[17]].n[17];

      /* delete upgoing velocities */
      lbfluid[k-next[4]].n[4]   = 0.0;
      lbfluid[k-next[10]].n[10] = 0.0;
      lbfluid[k-next[13]].n[13] = 0.0;
      lbfluid[k-next[14]].n[14] = 0.0;
      lbfluid[k-next[17]].n[17] = 0.0;

    }
    
  }

  /* top-down sweep */
  for (k=lblattice.halo_grid_volume-1-zperiod;k>=0;k--) {

    if (lbfluid[k].boundary == 1) {

      /* reflect downgoing velocities from higher indices */
      lbfluid[k].n[reflect[5]]  = lbfluid[k-next[5]].n[5];
      lbfluid[k].n[reflect[11]] = lbfluid[k-next[11]].n[11];
      lbfluid[k].n[reflect[12]] = lbfluid[k-next[12]].n[12];
      lbfluid[k].n[reflect[15]] = lbfluid[k-next[15]].n[15];
      lbfluid[k].n[reflect[16]] = lbfluid[k-next[16]].n[16];

      /* delete downgoing velocities */
      lbfluid[k-next[5]].n[5]   = 0.0;
      lbfluid[k-next[11]].n[11] = 0.0;
      lbfluid[k-next[12]].n[12] = 0.0;
      lbfluid[k-next[15]].n[15] = 0.0;
      lbfluid[k-next[16]].n[16] = 0.0;

    }
    
  }

}

/* canonic pre bounce */
MDINLINE void lb_slip_boundaries2() {

  int k;
  double *local_n, local_rho, local_j[3];
  double delta_j[3];

  double lb_slip_pref = lb_boundary_par.slip_pref;

  //int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   0;                       // ( 1, 0, 0) +
  next[1]  = - 0;                       // (-1, 0, 0)
  next[2]  =   0;                       // ( 0, 1, 0) +
  next[3]  = - 0;                       // ( 0,-1, 0)
  next[4]  = - 0;                       // ( 0, 0, 1)
  next[5]  =   0;                       // ( 0, 0,-1) +
  next[6]  =   0;                       // ( 1, 1, 0) +
  next[7]  = - 0;                       // (-1,-1, 0)
  next[8]  = - 0;                       // ( 1,-1, 0) 
  next[9]  =   0;                       // (-1, 1, 0) +
  next[10] = - 0;                       // ( 1, 0, 1)
  next[11] =   0;                       // (-1, 0,-1) +
  next[12] = - 0;                       // ( 1, 0,-1) +
  next[13] =   0;                       // (-1, 0, 1)
  next[14] = - 0;                       // ( 0, 1, 1)
  next[15] =   0;                       // ( 0,-1,-1) +
  next[16] = - 0;                       // ( 0, 1,-1) +
  next[17] =   0;                       // ( 0,-1, 1) 
  //int reflect[] = { 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 16, 17, 14, 15 };

  /* bottom-up sweep */
  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      //fprintf(stderr,"%d\n",k);

      local_n = lbfluid[k].n;
 
      /* calculate the density of upgoing velocities in the wall */
      local_rho =   local_n[4] 
                  + local_n[10] + local_n[13] 
	          + local_n[14] + local_n[17];

      /* calculate tangential component of upgoing velocities */
      local_j[0] =  local_n[10] - local_n[13];
      local_j[1] =  local_n[14] - local_n[17];

      //if (k==1714) {
      //	fprintf(stderr,"%d: (%.3e,%.3e,%.3e,%.3e,%.3e)\n",k,local_n[4],local_n[10],local_n[13],local_n[14],local_n[17]);
      //}

      //if (k==1597 || k==1714) {
      //	fprintf(stderr,"%d (%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f)\n",k,local_n[0],local_n[1],local_n[2],local_n[3],local_n[4],local_n[5],local_n[6],local_n[7],local_n[8],local_n[9],local_n[10],local_n[11],local_n[12],local_n[13],local_n[14],local_n[15],local_n[16],local_n[17]);
      //	fprintf(stderr,"%d (%.3f,%.3f,-)\n",k,local_j[0],local_j[1]);
      //}

      /* calculate friction at the wall (lattice units) */
      delta_j[0] = - lb_slip_pref * local_j[0]/local_rho * lbpar.tau;
      delta_j[1] = - lb_slip_pref * local_j[1]/local_rho * lbpar.tau;

      /* apply friction while reflecting back from the wall to lower indices */
      lbfluid[k-zperiod].n[0]  +=   1./6.*delta_j[0];
      lbfluid[k-zperiod].n[1]  += - 1./6.*delta_j[0];
      lbfluid[k-zperiod].n[2]  +=   1./6.*delta_j[1];
      lbfluid[k-zperiod].n[3]  += - 1./6.*delta_j[1];
      lbfluid[k-zperiod].n[4]  += 0.0;
      lbfluid[k-zperiod].n[5]  += local_n[4];
      lbfluid[k-zperiod].n[6]  +=   1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k-zperiod].n[7]  += - 1./12.*(delta_j[0]+delta_j[1]);
      lbfluid[k-zperiod].n[8]  +=   1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k-zperiod].n[9]  += - 1./12.*(delta_j[0]-delta_j[1]);
      lbfluid[k-zperiod].n[10] +=   1./12.*delta_j[0];
      lbfluid[k-zperiod].n[11] += local_n[13] - 1./12.*delta_j[0];
      lbfluid[k-zperiod].n[12] += local_n[10] + 1./12.*delta_j[0];
      lbfluid[k-zperiod].n[13] += - 1./12.*delta_j[0];
      lbfluid[k-zperiod].n[14] +=   1./12.*delta_j[1];
      lbfluid[k-zperiod].n[15] += local_n[17] - 1./12.*delta_j[1];
      lbfluid[k-zperiod].n[16] += local_n[14] + 1./12.*delta_j[1];
      lbfluid[k-zperiod].n[17] += - 1./12.*delta_j[1];

      //if (k==1714) {
      //	local_n = lbfluid[k-zperiod].n;
      //	fprintf(stderr,"%d: %e (%e,%e,-)\n",k,local_rho,local_j[0],local_j[1]);
      //	fprintf(stderr,"%d (%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)\n",k,local_n[0],local_n[1],local_n[2],local_n[3],local_n[4],local_n[5],local_n[6],local_n[7],local_n[8],local_n[9],local_n[10],local_n[11],local_n[12],local_n[13],local_n[14],local_n[15],local_n[16],local_n[17]);
      //}

      /* delete fluid in the wall */
      lbfluid[k].n[1] = 0.0;
      lbfluid[k].n[3] = 0.0;
      lbfluid[k].n[4] = 0.0;
      lbfluid[k].n[7] = 0.0;
      lbfluid[k].n[8] = 0.0;
      lbfluid[k].n[10] = 0.0;
      lbfluid[k].n[13] = 0.0;
      lbfluid[k].n[14] = 0.0;
      lbfluid[k].n[17] = 0.0;

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

}

/* canonic post bounce */
MDINLINE void lb_slip_boundaries3() {

  int k;
  double *local_n, local_rho, local_j[3];
  double delta_j[3];

  double lb_slip_pref = lb_boundary_par.slip_pref;

  //int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   0;                       // ( 1, 0, 0) +
  next[1]  = - 0;                       // (-1, 0, 0)
  next[2]  =   0;                       // ( 0, 1, 0) +
  next[3]  = - 0;                       // ( 0,-1, 0)
  next[4]  = - 0;                       // ( 0, 0, 1)
  next[5]  =   0;                       // ( 0, 0,-1) +
  next[6]  =   0;                       // ( 1, 1, 0) +
  next[7]  = - 0;                       // (-1,-1, 0)
  next[8]  = - 0;                       // ( 1,-1, 0) 
  next[9]  =   0;                       // (-1, 1, 0) +
  next[10] = - 0;                       // ( 1, 0, 1)
  next[11] =   0;                       // (-1, 0,-1) +
  next[12] = - 0;                       // ( 1, 0,-1) +
  next[13] =   0;                       // (-1, 0, 1)
  next[14] = - 0;                       // ( 0, 1, 1)
  next[15] =   0;                       // ( 0,-1,-1) +
  next[16] = - 0;                       // ( 0, 1,-1) +
  next[17] =   0;                       // ( 0,-1, 1) 
  //t reflect[] = { 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 16, 17, 14, 15 };

  /* bottom up sweep */
  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      //fprintf(stderr,"%d\n",k);

      local_n = lbfluid[k-zperiod].n;

      /* calculate the density of the wall reflected velocities */
      local_rho =   local_n[5] 
                  + local_n[11] + local_n[12] 
	          + local_n[15] + local_n[16];

      /* calculate tangential component of downgoing velocities */
      local_j[0] =  - local_n[11] + local_n[12];
      local_j[1] =  - local_n[15] + local_n[16];

      //if (k==1714) {
      //	fprintf(stderr,"%d: (%.3e,%.3e,%.3e,%.3e,%.3e)\n",k,local_n[5],local_n[11],local_n[12],local_n[15],local_n[16]);
      //}

      delta_j[0] = - lb_slip_pref * local_j[0]/local_rho * lbpar.tau;
      delta_j[1] = - lb_slip_pref * local_j[1]/local_rho * lbpar.tau;

      /* apply friction while reflecting back from the wall to lower indices */
      local_n[0]           +=   1./6.*delta_j[0];
      local_n[1]           += - 1./6.*delta_j[0];
      local_n[2]           +=   1./6.*delta_j[1];
      local_n[3]           += - 1./6.*delta_j[1];
      local_n[4]           +=   0.0;
      local_n[5]           +=   0.0;
      local_n[6]           +=   1./12.*(delta_j[0]+delta_j[1]);
      local_n[7]           += - 1./12.*(delta_j[0]+delta_j[1]);
      local_n[8]           +=   1./12.*(delta_j[0]-delta_j[1]);
      local_n[9]           += - 1./12.*(delta_j[0]-delta_j[1]);
      local_n[10]          +=   1./12.*delta_j[0];
      local_n[11]          += - 1./12.*delta_j[0];
      local_n[12]          +=   1./12.*delta_j[0];
      local_n[13]          += - 1./12.*delta_j[0];
      local_n[14]          +=   1./12.*delta_j[1];
      local_n[15]          += - 1./12.*delta_j[1];
      local_n[16]          +=   1./12.*delta_j[1];
      local_n[17]          += - 1./12.*delta_j[1];

      //if (k==1714) {
      //	fprintf(stderr,"%d: %e (%e,%e,-)\n",k,local_rho,local_j[0],local_j[1]);
      //	fprintf(stderr,"%d (%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)\n",k,local_n[0],local_n[1],local_n[2],local_n[3],local_n[4],local_n[5],local_n[6],local_n[7],local_n[8],local_n[9],local_n[10],local_n[11],local_n[12],local_n[13],local_n[14],local_n[15],local_n[16],local_n[17]);
      //}

    }
    
  }

  /* top down sweep */
  for (k=lblattice.halo_grid_volume-lblattice.halo_offset-1;k>=0;k--) {

    if (lbfluid[k].boundary) {

      //fprintf(stderr,"%d\n",k);

      local_n = lbfluid[k+zperiod].n;

      /* calculate the density of wall reflected velocities */
      local_rho =   local_n[4]
                  + local_n[10] + local_n[13]
	          + local_n[14] + local_n[17];

      /* calculate tangential component of upgoing */
      local_j[0] = local_n[10] - local_n[13];
      local_j[1] = local_n[14] - local_n[17];

      //if (k==1714) {
      //	fprintf(stderr,"%d: (%.3e,%.3e,%.3e,%.3e,%.3e)\n",k,local_n[4],local_n[10],local_n[13],local_n[14],local_n[17]);
      //}

      /* calculate friction at the wall (lattice units) */
      delta_j[0] = - lb_slip_pref * local_j[0]/local_rho * lbpar.tau;
      delta_j[1] = - lb_slip_pref * local_j[1]/local_rho * lbpar.tau;

      /* apply friction while reflecting back from the wall to higher indices */
      local_n[0]           +=   1./6.*delta_j[0];
      local_n[1]           += - 1./6.*delta_j[0];
      local_n[2]           +=   1./6.*delta_j[1];
      local_n[3]           += - 1./6.*delta_j[1];
      local_n[4]           += 0.0;
      local_n[5]           += 0.0;
      local_n[6]           +=   1./12.*(delta_j[0]+delta_j[1]);
      local_n[7]           += - 1./12.*(delta_j[0]+delta_j[1]);
      local_n[8]           +=   1./12.*(delta_j[0]-delta_j[1]);
      local_n[9]           += - 1./12.*(delta_j[0]-delta_j[1]);
      local_n[10]          +=   1./12.*delta_j[0];
      local_n[11]          += - 1./12.*delta_j[0];
      local_n[12]          +=   1./12.*delta_j[0];
      local_n[13]          += - 1./12.*delta_j[0];
      local_n[14]          +=   1./12.*delta_j[1];
      local_n[15]          += - 1./12.*delta_j[1];
      local_n[16]          +=   1./12.*delta_j[1];
      local_n[17]          += - 1./12.*delta_j[1];

      //if (k==1714) {
      //	fprintf(stderr,"%d: %e (%e,%e,-)\n",k,local_rho,local_j[0],local_j[1]);
      //	fprintf(stderr,"%d (%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)\n",k,local_n[0],local_n[1],local_n[2],local_n[3],local_n[4],local_n[5],local_n[6],local_n[7],local_n[8],local_n[9],local_n[10],local_n[11],local_n[12],local_n[13],local_n[14],local_n[15],local_n[16],local_n[17]);
      //}

    }
    
  }

}

/* simple post bounce */
MDINLINE void lb_slip_boundaries4() {

  int k;
  double *local_n, local_rho, local_j[3];
  double delta_j[3];

  double lb_slip_pref = lb_boundary_par.slip_pref;

  //int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   0;                       // ( 1, 0, 0) +
  next[1]  = - 0;                       // (-1, 0, 0)
  next[2]  =   0;                       // ( 0, 1, 0) +
  next[3]  = - 0;                       // ( 0,-1, 0)
  next[4]  = - 0;                       // ( 0, 0, 1)
  next[5]  =   0;                       // ( 0, 0,-1) +
  next[6]  =   0;                       // ( 1, 1, 0) +
  next[7]  = - 0;                       // (-1,-1, 0)
  next[8]  = - 0;                       // ( 1,-1, 0) 
  next[9]  =   0;                       // (-1, 1, 0) +
  next[10] = - 0;                       // ( 1, 0, 1)
  next[11] =   0;                       // (-1, 0,-1) +
  next[12] = - 0;                       // ( 1, 0,-1) +
  next[13] =   0;                       // (-1, 0, 1)
  next[14] = - 0;                       // ( 0, 1, 1)
  next[15] =   0;                       // ( 0,-1,-1) +
  next[16] = - 0;                       // ( 0, 1,-1) +
  next[17] =   0;                       // ( 0,-1, 1) 
  //t reflect[] = { 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 16, 17, 14, 15 };

  /* bottom up sweep */
  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume;k++) {

    if (lbfluid[k].boundary) {

      //fprintf(stderr,"%d\n",k);

      local_n = lbfluid[k-zperiod].n;

      /* calculate the density of the wall reflected velocities */
      local_rho =   local_n[5] 
                  + local_n[11] + local_n[12] 
	          + local_n[15] + local_n[16];

      /* calculate tangential component of downgoing velocities */
      local_j[0] =  - local_n[11] + local_n[12];
      local_j[1] =  - local_n[15] + local_n[16];

      delta_j[0] = - lb_slip_pref * local_j[0] * lbpar.tau;
      delta_j[1] = - lb_slip_pref * local_j[1] * lbpar.tau;

      /* apply friction while reflecting back from the wall to lower indices */
      local_n[11]          += - 0.5*delta_j[0];
      local_n[12]          +=   0.5*delta_j[0];
      local_n[15]          += - 0.5*delta_j[1];
      local_n[16]          +=   0.5*delta_j[1];

    }
    
  }

  /* top down sweep */
  for (k=lblattice.halo_grid_volume-lblattice.halo_offset-1;k>=0;k--) {

    if (lbfluid[k].boundary) {

      //fprintf(stderr,"%d\n",k);

      local_n = lbfluid[k+zperiod].n;

      /* calculate the density of wall reflected velocities */
      local_rho =   local_n[4]
                  + local_n[10] + local_n[13]
	          + local_n[14] + local_n[17];

      /* calculate tangential component of upgoing */
      local_j[0] = local_n[10] - local_n[13];
      local_j[1] = local_n[14] - local_n[17];

      //fprintf(stderr,"(%f,%f,-)\n",local_j[0],local_j[1]);

      /* calculate friction at the wall (lattice units) */
      delta_j[0] = - lb_slip_pref * local_j[0] * lbpar.tau;
      delta_j[1] = - lb_slip_pref * local_j[1] * lbpar.tau;

      /* apply friction while reflecting back from the wall to higher indices */
      local_n[10]          +=   0.5*delta_j[0];
      local_n[13]          += - 0.5*delta_j[0];
      local_n[14]          +=   0.5*delta_j[1];
      local_n[17]          += - 0.5*delta_j[1];

    }
    
  }

}

/* canonic without bounce */
MDINLINE void lb_slip_boundaries5() {

  int i,k;
  double *local_n, local_rho_up, local_rho_down, local_j_up[3], local_j_down[3];
  double delta_j_up[3], delta_j_down[3];
  double new_n[18];

  double lb_slip_pref = lb_boundary_par.slip_pref;

  //int yperiod = lblattice.halo_grid[0];
  //int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[18];
  next[0]  =   0;                       // ( 1, 0, 0) +
  next[1]  = - 0;                       // (-1, 0, 0)
  next[2]  =   0;                       // ( 0, 1, 0) +
  next[3]  = - 0;                       // ( 0,-1, 0)
  next[4]  = - 0;                       // ( 0, 0, 1)
  next[5]  =   0;                       // ( 0, 0,-1) +
  next[6]  =   0;                       // ( 1, 1, 0) +
  next[7]  = - 0;                       // (-1,-1, 0)
  next[8]  = - 0;                       // ( 1,-1, 0) 
  next[9]  =   0;                       // (-1, 1, 0) +
  next[10] = - 0;                       // ( 1, 0, 1)
  next[11] =   0;                       // (-1, 0,-1) +
  next[12] = - 0;                       // ( 1, 0,-1) +
  next[13] =   0;                       // (-1, 0, 1)
  next[14] = - 0;                       // ( 0, 1, 1)
  next[15] =   0;                       // ( 0,-1,-1) +
  next[16] = - 0;                       // ( 0, 1,-1) +
  next[17] =   0;                       // ( 0,-1, 1) 

  for (k=lblattice.halo_offset;k<lblattice.halo_grid_volume-lblattice.halo_offset;k++) {

    if (lbfluid[k].boundary) {

      local_n = lbfluid[k].n;
 
      /* calculate the density of upgoing velocities in the wall */
      local_rho_up   =   local_n[4] 
                       + local_n[10] + local_n[13] 
	               + local_n[14] + local_n[17];
      local_rho_down =   local_n[5]
                       + local_n[11] + local_n[12]
	               + local_n[15] + local_n[16];

      if (local_rho_up<ROUND_ERROR_PREC) local_rho_up=1.0;
      if (local_rho_down<ROUND_ERROR_PREC) local_rho_down=1.0;

      /* calculate tangential component of upgoing velocities */
      local_j_up[0]   =   local_n[10] - local_n[13];
      local_j_up[1]   =   local_n[14] - local_n[17];
      local_j_down[0] = - local_n[11] + local_n[12];
      local_j_down[1] = - local_n[15] + local_n[16];

      /* calculate friction at the wall (lattice units) */
      delta_j_up[0]   = - lb_slip_pref * local_j_up[0]/local_rho_up * lbpar.tau;
      delta_j_up[1]   = - lb_slip_pref * local_j_up[1]/local_rho_up * lbpar.tau;
      delta_j_down[0] = - lb_slip_pref * local_j_down[0]/local_rho_down * lbpar.tau;
      delta_j_down[1] = - lb_slip_pref * local_j_down[1]/local_rho_down * lbpar.tau;

      /* apply friction while reflecting back from the wall to lower indices */
      new_n[0]  = local_n[0] + 1./6.*delta_j_up[0] + 1./6.*delta_j_down[0];
      new_n[1]  = local_n[1] - 1./6.*delta_j_up[0] - 1./6.*delta_j_down[0];
      new_n[2]  = local_n[2] + 1./6.*delta_j_up[1] + 1./6.*delta_j_down[1];
      new_n[3]  = local_n[3] - 1./6.*delta_j_up[1] - 1./6.*delta_j_down[1];
      new_n[4]  = local_n[5];
      new_n[5]  = local_n[4];
      new_n[6]  = local_n[6] + 1./12.*(delta_j_up[0]+delta_j_up[1])
	                     + 1./12.*(delta_j_down[0]+delta_j_down[1]);
      new_n[7]  = local_n[7] - 1./12.*(delta_j_up[0]+delta_j_up[1])
                             - 1./12.*(delta_j_down[0]+delta_j_down[1]);
      new_n[8]  = local_n[8] + 1./12.*(delta_j_up[0]-delta_j_up[1])
                             + 1./12.*(delta_j_down[0]-delta_j_down[1]);
      new_n[9]  = local_n[9] - 1./12.*(delta_j_up[0]-delta_j_up[1])
                             - 1./12.*(delta_j_down[0]-delta_j_down[1]);
      new_n[10] = local_n[12] + 1./12.*delta_j_up[0] + 1./12.*delta_j_down[0];
      new_n[11] = local_n[13] - 1./12.*delta_j_up[0] - 1./12.*delta_j_down[0];
      new_n[12] = local_n[10] + 1./12.*delta_j_up[0] + 1./12.*delta_j_down[0];
      new_n[13] = local_n[11] - 1./12.*delta_j_up[0] - 1./12.*delta_j_down[0];
      new_n[14] = local_n[16] + 1./12.*delta_j_up[1] + 1./12.*delta_j_down[1];
      new_n[15] = local_n[17] - 1./12.*delta_j_up[1] - 1./12.*delta_j_down[1];
      new_n[16] = local_n[14] + 1./12.*delta_j_up[1] + 1./12.*delta_j_down[1];
      new_n[17] = local_n[15] - 1./12.*delta_j_up[1] - 1./12.*delta_j_down[1];

      //memcpy(local_n,new_n,18*sizeof(double));
      for (i=0;i<18;i++) {
	local_n[i] = new_n[i];
      }

    }
    
  }

}
#endif

MDINLINE void lb_boundary_friction(LB_FluidNode *node) {
    int i;

    double *local_n = node->n;
    double *local_j = node->j;
    double delta_j[3];

    double gamma = lb_boundary_par.slip_pref; 

    lb_calc_local_j(node);

    /* Einheiten? */
    
    delta_j[0] = - gamma * local_j[0] * lbpar.tau;
    delta_j[1] = - gamma * local_j[1] * lbpar.tau;
    delta_j[2] = 0.0;

#if 1
    local_n[1]  += 1./6.*delta_j[0];
    local_n[2]  -= 1./6.*delta_j[0];
    local_n[3]  += 1./6.*delta_j[1];
    local_n[4]  -= 1./6.*delta_j[1];
    local_n[5]  += 1./6.*delta_j[2];
    local_n[6]  -= 1./6.*delta_j[2];
    local_n[7]  += 1./12.*(delta_j[0]+delta_j[1]);
    local_n[8]  -= 1./12.*(delta_j[0]+delta_j[1]);
    local_n[9]  += 1./12.*(delta_j[0]-delta_j[1]);
    local_n[10] -= 1./12.*(delta_j[0]-delta_j[1]);
    local_n[11] += 1./12.*(delta_j[0]+delta_j[2]);
    local_n[12] -= 1./12.*(delta_j[0]+delta_j[2]);
    local_n[13] += 1./12.*(delta_j[0]-delta_j[2]);
    local_n[14] -= 1./12.*(delta_j[0]-delta_j[2]);
    local_n[15] += 1./12.*(delta_j[1]+delta_j[2]);
    local_n[16] -= 1./12.*(delta_j[1]+delta_j[2]);
    local_n[17] += 1./12.*(delta_j[1]-delta_j[2]);
    local_n[18] -= 1./12.*(delta_j[1]-delta_j[2]);  
#else
    local_n[11]  += 1./1.*delta_j[0];
    //local_n[12]  -= 1./4.*delta_j[0];
    //local_n[13]  += 1./4.*delta_j[0];
    //local_n[14]  -= 1./4.*delta_j[0];
    local_n[15]  += 1./1.*delta_j[1];
    //local_n[16]  -= 1./4.*delta_j[1];
    //local_n[17]  += 1./4.*delta_j[1];
    //local_n[18]  -= 1./4.*delta_j[1];
#endif
    for (i=0;i<19;i++) {
	if (local_n[i]<0.0) {
	    fprintf(stderr,"%d: negative population %d\n",this_node,i);
	}
    }

}

MDINLINE void lb_boundary_equilibrium1(LB_FluidNode *node, const double local_rho, const double *local_j, const double *local_pi) {
  int i, j, k, l, m;

  /* model dependent variables */
  double (*c)[3]       = lbmodel.c;
  double (*lbcoeff)[4] = lbmodel.coeff;
  double *w            = lbmodel.w;
  double c_sound_sq    = lbmodel.c_sound_sq;

  /* boundary dependent variables */
  //int mask     = node->boundary;
  double *nvec = node->nvec;

  /* local variables */
  double *local_n = node->n;
  double trace, tmp;
  double uvec[3], vvec[3];
  double j_par[3], j_perp, j_abs;
  double pi_uu, pi_vv, pi_uv, pi_un, pi_vn;
  double coeff[7];
  double *lambda;
  double *b;
  double **A;
  int *perms;

  /* calculate the actual hydrodynamic fields */
  //lb_calc_local_fields(node,1);

  for (i=0; i<lbmodel.n_veloc; i++) {
    //if(scalar(c[i],nvec) < 0.0)
    {
      fprintf(stderr,"local_n[%d]=%.24e\n",i,local_n[i]+lbmodel.coeff[i][0]*lbpar.rho_lb_units);
    }
  }

  fprintf(stderr,"Start: local_rho=%f local_j=(%e,%e,%e) Pi=(%.12e,%e,%.12e,%e,%f,%.12e)\n",local_rho,local_j[0],local_j[1],local_j[2],local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]);
  
  //if (node->boundary == 1) {
  //  local_n[5] += local_n[6];
  //  local_n[11] += local_n[13];
  //  local_n[14] += local_n[12];
  //  local_n[15] += local_n[17];
  //  local_n[18] += local_n[16];
  //  local_n[6] = 0.;
  //  local_n[13] = 0.;
  //  local_n[12] = 0.;
  //  local_n[17] = 0.;
  //  local_n[16] = 0.;
  //}
  //if (node->boundary == -1) {
  //  local_n[6] += local_n[5];
  //  local_n[13] += local_n[11];
  //  local_n[12] += local_n[14];
  //  local_n[17] += local_n[15];
  //  local_n[16] += local_n[18];
  //  local_n[5] = 0.;
  //  local_n[11] = 0.;
  //  local_n[14] = 0.;
  //  local_n[15] = 0.;
  //  local_n[18] = 0.;
  //}
  //
  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  //if(scalar(c[i],nvec) < 0.0)
  //  {
  //    fprintf(stderr,"post local_n[%d]=%.24e\n",i,local_n[i]);
  //  }
  //}
  //
  //lb_calc_local_fields(node,1);
  //
  //fprintf(stderr,"Start: local_rho=%f local_j=(%e,%e,%e) Pi=(%.12e,%e,%.12e,%e,%f,%.12e)\n",local_rho,local_j[0],local_j[1],local_j[2],local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]);
  //
  //return;
  
  double rho = local_rho;
  double avg_rho = lbpar.rho_lb_units;

  double v[3] = { 0., 0., 0. };
  v[0] = local_j[0];
  v[1] = local_j[1];
  v[2] = 0.0;

  /* \TODO: Einheiten und fuer beliebige Orientierung */
  double gamma = lb_boundary_par.slip_pref;
  double tau = lbpar.tau;
  double delta_j[3];
  delta_j[0] = - gamma * v[0] * tau;
  delta_j[1] = - gamma * v[1] * tau;
  delta_j[2] = - gamma * v[2] * tau;

  v[0] += delta_j[0];
  v[1] += delta_j[1];
  v[2] += delta_j[2];

  double pi[6] = { 0., 0., 0., 0., 0., 0. };
  //pi[0] = rho*c_sound_sq+rho*v[0]*v[0];
  //pi[1] = rho*v[0]*v[1];
  //pi[2] = rho*c_sound_sq+rho*v[1]*v[1];
  //pi[3] = rho*v[0]*v[2];
  //pi[4] = rho*v[1]*v[2];
  //pi[5] = rho*c_sound_sq+rho*v[2]*v[2];

  pi[0] = local_pi[0]; //-trace/2.0;
  pi[1] = local_pi[1];
  pi[2] = local_pi[2]; //-trace/2.0;
  pi[3] = local_pi[3];// + node->boundary*delta_j[0];
  pi[4] = local_pi[4];// + node->boundary*delta_j[1];  
  pi[5] = 2.0*local_pi[5];

  fprintf(stderr,"use for bulk: rho=%f v=(%.12e,%f,%f) ",rho,v[0],v[1],v[2]);
  fprintf(stderr,"pi=(%f,%e,%f,%e,%f,%e)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);

  double rhoc_sq = local_rho*lbmodel.c_sound_sq;
  pi[0] -= rhoc_sq;
  pi[2] -= rhoc_sq;
  pi[5] -= rhoc_sq;

  trace = pi[0] + pi[2] + pi[5];

  /* loop over all links to do the bulk-like part */
  for (i=0; i<lbmodel.n_veloc; i++) {

    if (scalar(c[i],nvec) > 0.0) 
    {
      /* fluid pointing links are treated like in the bulk */

      //local_n[i] = rho*( 
      //			lbcoeff[i][0] 
      //			+ lbcoeff[i][1] * scalar(v,c[i])
      //			+ lbcoeff[i][2] * scalar(v,c[i])*scalar(v,c[i])
      //			+ lbcoeff[i][3] * scalar(v,v)
      //		 );

      //int reflect[] = { 0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 13, 14, 11, 12, 17, 18, 15, 16 }; 
      //local_n[i] = local_n[reflect[i]];

      local_n[i]  = lbcoeff[i][0] * (rho - avg_rho);
      //fprintf(stderr,"local_n[%d]=%.12e rho-avg_rho=%e\n",i,local_n[i],rho-avg_rho);
      local_n[i] += lbcoeff[i][1] * scalar(v,c[i]);
      //fprintf(stderr,"local_n[%d]=%.12e scalar=%e\n",i,local_n[i],lbcoeff[i][1]*scalar(v,c[i]));

      tmp = 0.0;
      m = 0;
      for (k=0;k<3;k++) {
      	for (l=0;l<k;l++) {
      	  /* non-diagonal elements */
      	  tmp += 2.0 * pi[m] * c[i][k] * c[i][l];
      	  m++;
      	}
      	/* diagonal elements */
      	tmp += pi[m] * c[i][k] * c[i][k];
      	m++;
      }
      
      local_n[i] += lbcoeff[i][2] * tmp;
      
      local_n[i] += lbcoeff[i][3] * trace;

      fprintf(stderr,"local_n[%d]=%.12e tmp=%e trace=%e\n",i,local_n[i]+lbcoeff[i][0]*avg_rho,tmp,trace);
    } else {
      /* all other links first get zero */
      local_n[i] = -lbcoeff[i][0]*avg_rho;
    }

  }

  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  if (scalar(c[i],nvec) <= 0.0) {
  //    local_n[i] = 0.0;
  //  }
  //}

  rho = avg_rho;
  v[0] = v[1] = v[2] = 0.0;
  pi[0] = pi[2] = pi[5] = avg_rho/3.0;
  pi[1] = pi[3] = pi[4] = 0.0;

  for (i=0; i<lbmodel.n_veloc; i++) {
      tmp = local_n[i];
      rho   += tmp;
      v[0]  += tmp*c[i][0];
      v[1]  += tmp*c[i][1];
      v[2]  += tmp*c[i][2];
      //fprintf(stderr,"v=(%e,%e,%e)\n",v[0],v[1],v[2]);
      pi[0] += tmp*c[i][0]*c[i][0];
      pi[1] += tmp*c[i][0]*c[i][1];
      pi[2] += tmp*c[i][1]*c[i][1];
      pi[3] += tmp*c[i][0]*c[i][2];
      pi[4] += tmp*c[i][1]*c[i][2];
      pi[5] += tmp*c[i][2]*c[i][2];
  }

  //v[0] /= rho;
  //v[1] /= rho;
  //v[2] /= rho;

  fprintf(stderr,"after bulk: rho=%f\tj=(%.12e,%f,%.12e)\t",rho,v[0],v[1],v[2]);
  fprintf(stderr,"pi=(%f,%f,%f,%e,%f,%f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);

  /* Einheiten und rho! */

  rho = local_rho - rho;

  v[0] = local_j[0] - v[0];
  v[1] = local_j[1] - v[1];
  v[2] = local_j[2] - v[2];

  pi[0] = local_pi[0] - pi[0];
  pi[1] = local_pi[1] - pi[1];
  pi[2] = local_pi[2] - pi[2];
  pi[3] = local_pi[3] - pi[3];
  pi[4] = local_pi[4] - pi[4];
  pi[5] = local_pi[5] - pi[5];

  v[0] += delta_j[0];
  v[1] += delta_j[1];
  v[2] += delta_j[2];

  fprintf(stderr,"diff: rho=%f\tj=(%e,%e,%e)\t",rho,v[0],v[1],v[2]);
  fprintf(stderr,"pi=(%f,%e,%f,%e,%f,%f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]);

  /* the difference has to be assigned to the surface links */

  /* calculate basis vectors from boundary orientation and flow velocity */
  j_perp = scalar(v,nvec);
  j_par[0] = v[0] - j_perp*nvec[0];
  j_par[1] = v[1] - j_perp*nvec[1];
  j_par[2] = v[2] - j_perp*nvec[2];

  j_abs = sqrt(sqrlen(j_par));
  if (j_abs) {
    unit_vector(j_par,uvec);
  } else {
    uvec[0] = 0.0;
    uvec[1] = -nvec[2];
    uvec[2] =  nvec[1];
  }
  
  vector_product(nvec,uvec,vvec);

  pi_uu = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      /* nondiagonal components */
      pi_uu += 2*pi[k]*uvec[i]*uvec[j];
      k++;
    }
    /* diagonal components */
    pi_uu += pi[k]*uvec[i]*uvec[i];
    k++;
  }
  
  pi_vv = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      /* nondiagonal components */
      pi_vv += 2*pi[k]*vvec[i]*vvec[j];
      k++;
    }
    /* diagonal components */
    pi_vv += pi[k]*vvec[i]*vvec[i];
    k++;
  }

  pi_uv = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      /* nondiagonal components */
      pi_uv += pi[k]*uvec[i]*vvec[j] + pi[k]*uvec[j]*vvec[i];
      k++;
    }
    pi_uv += pi[k]*uvec[i]*vvec[i];
    k++;
  }

  pi_un = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      pi_un += pi[k]*uvec[i]*nvec[j] + pi[k]*uvec[j]*nvec[i];
      k++;
    }
    pi_un += pi[k]*uvec[i]*nvec[i];
    k++;
  }

  pi_vn = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      pi_vn += pi[k]*vvec[i]*nvec[j] + pi[k]*vvec[j]*nvec[i];
      k++;
    }
    pi_vn += pi[k]*vvec[i]*nvec[i];
    k++;
  }

  //fprintf(stderr,"n=(%e,%e,%e)\n",nvec[0],nvec[1],nvec[2]);
  //fprintf(stderr,"u=(%e,%e,%e)\n",uvec[0],uvec[1],uvec[2]);
  //fprintf(stderr,"v=(%e,%e,%e)\n",vvec[0],vvec[1],vvec[2]);
  //fprintf(stderr,"j_abs=%.12e\t",j_abs);
  //fprintf(stderr,"pi_uu=%e pi_vv=%e pi_uv=%e\n",pi_uu,pi_vv,pi_uv);

  int n_constraints = 6;

  /* memory allocation */
  perms = malloc(n_constraints*sizeof(int));
  b = malloc(n_constraints*sizeof(double));
  A = malloc(n_constraints*sizeof(double *));
  *A = malloc(n_constraints*n_constraints*sizeof(double));

  /* store the desired values of the moments in a vector */
  b[0] = rho;
  b[1] = j_abs;
  b[2] = 0.0;
  b[3] = pi_uu;//(local_rho) ? local_rho*c_sound_sq + SQR(j_abs)/(local_rho) : 0.0;
  b[4] = pi_vv;//local_rho*c_sound_sq;
  b[5] = pi_uv;//0.0;
  //b[6] = pi_un;
  //b[7] = pi_vn;
  
  /* initialize matrix */
  for (i=0; i<n_constraints; i++) {
    A[i] = A[0]+i*n_constraints;
    for (j=0; j<n_constraints; j++) {
      A[i][j] = 0.0;
    }
  }

  /* loop over all links */
  for (i=0; i<lbmodel.n_veloc; i++) {

    /* only boundary surface links contribute */
    if (scalar(c[i],nvec) == 0.0) {
      //if (!(mask & (1 << i))) {
      //fprintf(stderr,"%d\n",i);
      /* calculate the coefficients of the Lagrange multipliers */
      coeff[0] = 1.0 ;
      coeff[1] = scalar(c[i],uvec);
      //coeff[2] = scalar(c[i],nvec);
      coeff[2] = scalar(c[i],vvec);
      coeff[3] = SQR(scalar(c[i],uvec));
      coeff[4] = SQR(scalar(c[i],vvec));
      coeff[5] = scalar(c[i],uvec)*scalar(c[i],vvec);
      //coeff[6] = scalar(c[i],uvec)*scalar(c[i],nvec);
      //coeff[7] = scalar(c[i],vvec)*scalar(c[i],nvec);

      /* sum up the contributions for the moments */
      for (j=0; j<n_constraints; j++) {
	for (k=0; k<n_constraints; k++) {
	  A[k][j] += w[i]*coeff[j]*coeff[k];
	  //A[1][j] += w[i]*coeff[j]*coeff[1];
	  //A[2][j] += w[i]*coeff[j]*coeff[2];
	  //A[3][j] += w[i]*coeff[j]*coeff[3];
	  //A[4][j] += w[i]*coeff[j]*coeff[4];
	  //A[5][j] += w[i]*coeff[j]*coeff[5];
	  //A[6][j] += w[i]*coeff[j]*coeff[6];
	}
      }

    }

  }
  
  for (i=0; i<n_constraints; i++) {
    for (j=0; j<n_constraints; j++) {
      fprintf(stderr,"%f ",A[i][j]);
    }
    fprintf(stderr,"\t%.12e\n",b[i]);
  }

  /* solve the linear equation system */
  if (lu_decompose_matrix(A,n_constraints,perms)) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{110 singular matrix while solving constraints in lb_boundary_equilibrium()} ");
    fprintf(stderr, "Matrix singular!\n");
    return;
  }
  lu_solve_system(A,n_constraints,perms,b);
  
  lambda = b;

  //fprintf(stderr,"lambda=(%f,%.12e,%f,%f,%f,%f)\n",lambda[0],lambda[1],lambda[2],lambda[3],lambda[4],lambda[5]);

  /* loop over all links */
  for (i=0; i<lbmodel.n_veloc; i++) {

    /* only fluid links have to be calculated */
    //if (mask & (1 << i)) {
    //  local_n[i] = 0.0;
    //}
    //else {
    if (scalar(c[i],nvec) == 0.0) {

      /* calculate the populations of the fluid links */
      local_n[i] = w[i] * ( lambda[0]
			    + lambda[1] * scalar(c[i],uvec)
			    //+ lambda[2] * scalar(c[i],nvec)
			    + lambda[2] * scalar(c[i],vvec)
			    + lambda[3] * SQR(scalar(c[i],uvec))
			    + lambda[4] * SQR(scalar(c[i],vvec))
			    + lambda[5] * scalar(c[i],uvec)*scalar(c[i],vvec)
			    //+ lambda[6] * scalar(c[i],uvec)*scalar(c[i],nvec)
			    //+ lambda[7] * scalar(c[i],vvec)*scalar(c[i],nvec)
			    );

      //fprintf(stderr,"%f * (%f,%f,%f,%f,%f,%f)\n",w[i],lambda[0],lambda[1]*scalar(c[i],uvec),lambda[2]*scalar(c[i],vvec),lambda[3]*SQR(scalar(c[i],uvec)),lambda[4]*SQR(scalar(c[i],vvec)),lambda[5]*scalar(c[i],uvec)*scalar(c[i],vvec));
      //fprintf(stderr,"local_n[%d]=%f\n",i,local_n[i]);

      if (local_n[i] < 0.0) {
	fprintf(stderr,"Negative population local_n[%d]=%f\n",i,local_n[i]);
      }
      if (isnan(local_n[i])) {
	fprintf(stderr,"Population not a number local_n[%d]=%f\n",i,local_n[i]);
	errexit();
      }
      local_n[i] -= lbcoeff[i][0]*avg_rho;
      fprintf(stderr,"local_n[%d]=%e\n",i,local_n[i]+lbcoeff[i][0]*avg_rho);      
    } else {
      //local_n[i] = lbcoeff[i][0]*avg_rho;
      fprintf(stderr,"local_n[%d]=%e\n",i,local_n[i]+lbcoeff[i][0]*avg_rho);
    }

  }



  lb_calc_local_fields(node,1);
  fprintf(stderr,"Final: rho=%.12e\t",*(node->rho));
  fprintf(stderr,"j=(%.12e,%f,%.12e)\t",node->j[0],node->j[1],node->j[2]);
  fprintf(stderr,"pi=(%f,%f,%f,%e,%f,%f)\n",node->pi[0],node->pi[1],node->pi[2],node->pi[3],node->pi[4],node->pi[5]);

  //if (fabs(node->pi[3])>ROUND_ERROR_PREC) {
  //  fprintf(stderr,"x-z-shear %e!\n",node->pi[3]);
  //  errexit();
  //}

  /* free memory */
  free(*A);
  free(A);
  free(b);
  free(perms);


}


extern double lblambda;
extern double lblambda_bulk;

MDINLINE void lb_boundary_collisions(LB_FluidNode *local_node) {
  int i;

  const double local_rho = *(local_node->rho);
  double *local_j  = local_node->j;
  double *local_pi = local_node->pi;
  double local_pi_eq[6];
  double trace, trace_eq;
  double tmp;

  const double rhoc_sq = local_rho*lbmodel.c_sound_sq;
  const double onepluslambda = 1.0 + lblambda;

  /* specular reflections */
  local_j[2] = - local_j[2];
  
  local_pi[3] = - local_pi[3];
  local_pi[4] = - local_pi[4];

  lb_boundary_equilibrium2(local_node, local_rho, local_j, local_pi);

  double save_pi[6];
  for (i=0;i<6;i++) save_pi[i] = local_pi[i];

  /* calculate the equilibrium pressure tensor */
  local_pi_eq[0] = rhoc_sq + local_j[0]*local_j[0]/local_rho;
  tmp = local_j[1]/local_rho;
  local_pi_eq[1] = local_j[0]*tmp;
  local_pi_eq[2] = rhoc_sq + local_j[1]*tmp;
  tmp = local_j[2]/local_rho;
  local_pi_eq[3] = local_j[0]*tmp;
  local_pi_eq[4] = local_j[1]*tmp;
  local_pi_eq[5] = local_node->boundary*local_j[2];

  lb_boundary_equilibrium2(local_node, local_rho, local_j, local_pi_eq);

  for (i=0;i<6;i++) local_pi[i] = save_pi[i]; 

  fprintf(stderr,"onepluslambda=%f\tlambda_v=%f\n",onepluslambda,lblambda_bulk);
  fprintf(stderr,"pi=(%f,%f,%f,%f,%f,%f)\n",local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]);

  fprintf(stderr,"pi_eq=(%f,%f,%f,%f,%f,%f)\n",local_pi_eq[0],local_pi_eq[1],local_pi_eq[2],local_pi_eq[3],local_pi_eq[4],local_pi_eq[5]);

  /* calculate the traces */
  trace_eq = local_pi_eq[0] + local_pi_eq[2] + local_pi_eq[5];
  trace = local_pi[0] + local_pi[2] + local_pi[5];
    
  /* relax the local pressure tensor */
  local_pi[0] = local_pi_eq[0] + onepluslambda*(local_pi[0] - local_pi_eq[0]);
  local_pi[1] = local_pi_eq[1] + onepluslambda*(local_pi[1] - local_pi_eq[1]);
  local_pi[2] = local_pi_eq[2] + onepluslambda*(local_pi[2] - local_pi_eq[2]);
  local_pi[3] = local_pi_eq[3] + onepluslambda*(local_pi[3] - local_pi_eq[3]);
  local_pi[4] = local_pi_eq[4] + onepluslambda*(local_pi[4] - local_pi_eq[4]);
  local_pi[5] = local_pi_eq[5] + onepluslambda*(local_pi[5] - local_pi_eq[5]);  
  tmp = 1./3.*(lblambda_bulk-lblambda)*(trace - trace_eq);
  local_pi[0] += tmp;
  local_pi[2] += tmp;
  local_pi[5] += tmp;

  fprintf(stderr,"pi*=(%f,%f,%f,%f,%f,%f)\n",local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]);
  
}

MDINLINE void lb_boundary_equilibrium2(LB_FluidNode *node, const double local_rho, const double *local_j, const double *local_pi) {
  int i, j, k, l, m;

  /* model dependent variables */
  double (*c)[3]       = lbmodel.c;
  double (*lbcoeff)[4] = lbmodel.coeff;
  double *w            = lbmodel.w;
  double c_sound_sq    = lbmodel.c_sound_sq;

  /* boundary dependent variables */
  //int mask     = node->boundary;
  double *nvec = node->nvec;

  /* local variables */
  double *local_n = node->n;
  double trace, tmp;
  double uvec[3], vvec[3];
  double j_par[3], j_perp, j_abs;
  double pi_uu, pi_vv, pi_uv, pi_un, pi_vn;
  double coeff[7];
  double *lambda;
  double *b;
  double **A;
  int *perms;

  /* calculate the actual hydrodynamic fields */
  //lb_calc_local_fields(node,1);

  //for (i=0; i<lbmodel.n_veloc; i++) {
  //  //if(scalar(c[i],nvec) < 0.0)
  //  {
  //    fprintf(stderr,"local_n[%d]=%.24e\n",i,local_n[i]);
  //  }
  //}

  LB_TRACE(fprintf(stderr,"Start: local_rho=%f local_j=(%e,%e,%e) Pi=(%.12e,%e,%.12e,%e,%f,%.12e)\n",local_rho,local_j[0],local_j[1],local_j[2],local_pi[0],local_pi[1],local_pi[2],local_pi[3],local_pi[4],local_pi[5]));
  
  double rho = local_rho;
  double avg_rho = lbpar.rho_lb_units;

  double v[3] = { 0., 0., 0. };
  v[0] = local_j[0];
  v[1] = local_j[1];
  v[2] = -local_j[2];

  v[2] = node->boundary*1./6.;

  /* \TODO: Einheiten und fuer beliebige Orientierung */
  double gamma = lb_boundary_par.slip_pref;
  double tau = lbpar.tau;
  double delta_j[3];
  delta_j[0] = - gamma * v[0] * tau;
  delta_j[1] = - gamma * v[1] * tau;
  delta_j[2] = - gamma * v[2] * tau;

  double pi[6] = { 0., 0., 0., 0., 0., 0. };
  //pi[0] = rho*c_sound_sq+rho*v[0]*v[0];
  //pi[1] = rho*v[0]*v[1];
  //pi[2] = rho*c_sound_sq+rho*v[1]*v[1];
  //pi[3] = rho*v[0]*v[2];
  //pi[4] = rho*v[1]*v[2];
  //pi[5] = rho*c_sound_sq+rho*v[2]*v[2];

  pi[0] = local_pi[0];// + 2.*v[0]*delta_j[0] + SQR(delta_j[0]);
  pi[1] = local_pi[1];// + v[0]*delta_j[1] + v[1]*delta_j[0] + delta_j[0]*delta_j[1];
  pi[2] = local_pi[2];// + 2.*v[1]*delta_j[1] + SQR(delta_j[1]);
  pi[3] = local_pi[3];// + v[0]*delta_j[2] + v[2]*delta_j[0] + delta_j[0]*delta_j[2];
  pi[4] = local_pi[4];// + v[1]*delta_j[2] + v[2]*delta_j[1] + delta_j[1]*delta_j[2];
  pi[5] = local_pi[5];// + 2.*v[2]*delta_j[2] + SQR(delta_j[2]);

  //pi[3] =   local_n[11] + local_n[12] - local_n[13] - local_n[14];
  //pi[4] =   local_n[15] + local_n[16] - local_n[17] - local_n[18];
  //
  //pi[3] -= node->boundary*delta_j[0];
  //pi[4] -= node->boundary*delta_j[1];

  //v[0] += delta_j[0];
  //v[1] += delta_j[1];
  //v[2] += delta_j[2];

  LB_TRACE(fprintf(stderr,"use: rho=%f\tj=(%e,%e,%e)\t",rho,v[0],v[1],v[2]));
  LB_TRACE(fprintf(stderr,"pi=(%f,%e,%f,%e,%f,%f)\n",pi[0],pi[1],pi[2],pi[3],pi[4],pi[5]));

  /* calculate basis vectors from boundary orientation and flow velocity */
  j_perp = scalar(v,nvec);
  j_par[0] = v[0] - j_perp*nvec[0];
  j_par[1] = v[1] - j_perp*nvec[1];
  j_par[2] = v[2] - j_perp*nvec[2];

  j_abs = sqrt(sqrlen(j_par));
  if (j_abs) {
    unit_vector(j_par,uvec);
  } else {
    uvec[0] = 0.0;
    uvec[1] = -nvec[2];
    uvec[2] =  nvec[1];
  }
  
  vector_product(nvec,uvec,vvec);

  pi_uu = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      /* nondiagonal components */
      pi_uu += 2*pi[k]*uvec[i]*uvec[j];
      k++;
    }
    /* diagonal components */
    pi_uu += pi[k]*uvec[i]*uvec[i];
    k++;
  }
  
  pi_vv = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      /* nondiagonal components */
      pi_vv += 2*pi[k]*vvec[i]*vvec[j];
      k++;
    }
    /* diagonal components */
    pi_vv += pi[k]*vvec[i]*vvec[i];
    k++;
  }

  pi_uv = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      /* nondiagonal components */
      pi_uv += pi[k]*uvec[i]*vvec[j] + pi[k]*uvec[j]*vvec[i];
      k++;
    }
    pi_uv += pi[k]*uvec[i]*vvec[i];
    k++;
  }

  pi_un = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      pi_un += pi[k]*uvec[i]*nvec[j] + pi[k]*uvec[j]*nvec[i];
      k++;
    }
    pi_un += pi[k]*uvec[i]*nvec[i];
    k++;
  }

  pi_vn = 0.0;
  k = 0;
  for (i=0; i<3; i++) {
    for (j=0; j<i; j++) {
      pi_vn += pi[k]*vvec[i]*nvec[j] + pi[k]*vvec[j]*nvec[i];
      k++;
    }
    pi_vn += pi[k]*vvec[i]*nvec[i];
    k++;
  }

  //j_perp = 1./6.;//local_j[2];

  //fprintf(stderr,"n=(%e,%e,%e)\n",nvec[0],nvec[1],nvec[2]);
  //fprintf(stderr,"u=(%e,%e,%e)\n",uvec[0],uvec[1],uvec[2]);
  //fprintf(stderr,"v=(%e,%e,%e)\n",vvec[0],vvec[1],vvec[2]);
  //fprintf(stderr,"j_abs=%.12e ",j_abs);
  //fprintf(stderr,"j_perp=.12e ",j_perp);
  //fprintf(stderr,"pi_uu=%e pi_vv=%e pi_uv=%e pi_un=%e pi_vn=%e\n",pi_uu,pi_vv,pi_uv,pi_un,pi_vn);

  int n_constraints = 7;

  /* memory allocation */
  perms = malloc(n_constraints*sizeof(int));
  b = malloc(n_constraints*sizeof(double));
  A = malloc(n_constraints*sizeof(double *));
  *A = malloc(n_constraints*n_constraints*sizeof(double));

  /* store the desired values of the moments in a vector */
  b[0] = rho;
  b[1] = j_abs;
  b[2] = 0.0;
  b[3] = j_perp;
  b[4] = pi_uu;//(local_rho) ? local_rho*c_sound_sq + SQR(j_abs)/(local_rho) : 0.0;
  b[5] = pi_vv;//local_rho*c_sound_sq;
  b[6] = pi_uv;//0.0;
  //b[7] = pi_un;
  //b[8] = pi_vn;
  
  /* initialize matrix */
  for (i=0; i<n_constraints; i++) {
    A[i] = A[0]+i*n_constraints;
    for (j=0; j<n_constraints; j++) {
      A[i][j] = 0.0;
    }
  }

  /* loop over all links */
  for (i=0; i<lbmodel.n_veloc; i++) {

    /* only boundary surface links contribute */
    if (scalar(c[i],nvec) >= 0.0) {
      //if (!(mask & (1 << i))) {
      //fprintf(stderr,"%d\n",i);
      /* calculate the coefficients of the Lagrange multipliers */
      coeff[0] = 1.0 ;
      coeff[1] = scalar(c[i],uvec);
      coeff[2] = scalar(c[i],vvec);
      coeff[3] = scalar(c[i],nvec);
      coeff[4] = SQR(scalar(c[i],uvec));
      coeff[5] = SQR(scalar(c[i],vvec));
      coeff[6] = scalar(c[i],uvec)*scalar(c[i],vvec);
      //coeff[7] = scalar(c[i],uvec)*scalar(c[i],nvec);
      //coeff[8] = scalar(c[i],vvec)*scalar(c[i],nvec);

      /* sum up the contributions for the moments */
      for (j=0; j<n_constraints; j++) {
	for (k=0; k<n_constraints; k++) {
	  A[k][j] += w[i]*coeff[j]*coeff[k];
	  //A[1][j] += w[i]*coeff[j]*coeff[1];
	  //A[2][j] += w[i]*coeff[j]*coeff[2];
	  //A[3][j] += w[i]*coeff[j]*coeff[3];
	  //A[4][j] += w[i]*coeff[j]*coeff[4];
	  //A[5][j] += w[i]*coeff[j]*coeff[5];
	  //A[6][j] += w[i]*coeff[j]*coeff[6];
	}
      }

    }

  }
  
  //for (i=0; i<n_constraints; i++) {
  //  for (j=0; j<n_constraints; j++) {
  //    fprintf(stderr,"%f ",A[i][j]);
  //  }
  //  fprintf(stderr,"\t%.12e\n",b[i]);
  //}

  /* solve the linear equation system */
  if (lu_decompose_matrix(A,n_constraints,perms)) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{110 singular matrix while solving constraints in lb_boundary_equilibrium()} ");
    fprintf(stderr, "Matrix singular!\n");
    return;
  }
  lu_solve_system(A,n_constraints,perms,b);
  
  lambda = b;

  //fprintf(stderr,"lambda=(%e,%e,%e,%e,%e,%e,%e,%e,%e)\n",lambda[0],lambda[1],lambda[2],lambda[3],lambda[4],lambda[5],lambda[6],lambda[7],lambda[8]);

  /* loop over all links */
  for (i=0; i<lbmodel.n_veloc; i++) {

    /* only fluid links have to be calculated */
    //if (mask & (1 << i)) {
    //  local_n[i] = 0.0;
    //}
    //else {
    if (scalar(c[i],nvec) >= 0.0) {

      /* calculate the populations of the fluid links */
      local_n[i] = w[i] * ( lambda[0]
			    + lambda[1] * scalar(c[i],uvec)
			    + lambda[2] * scalar(c[i],vvec)
			    + lambda[3] * scalar(c[i],nvec)
			    + lambda[4] * SQR(scalar(c[i],uvec))
			    + lambda[5] * SQR(scalar(c[i],vvec))
			    + lambda[6] * scalar(c[i],uvec)*scalar(c[i],vvec)
			    //+ lambda[7] * scalar(c[i],uvec)*scalar(c[i],nvec)
			    //+ lambda[8] * scalar(c[i],vvec)*scalar(c[i],nvec)
			    );

      //fprintf(stderr,"%f * (%f,%f,%f,%f,%f,%f)\n",w[i],lambda[0],lambda[1]*scalar(c[i],uvec),lambda[2]*scalar(c[i],vvec),lambda[3]*SQR(scalar(c[i],uvec)),lambda[4]*SQR(scalar(c[i],vvec)),lambda[5]*scalar(c[i],uvec)*scalar(c[i],vvec));
      //fprintf(stderr,"local_n[%d]=%f\n",i,local_n[i]);

      //if (local_n[i] < 0.0) {
      //	fprintf(stderr,"Negative population local_n[%d]=%f\n",i,local_n[i]);
      //}
      if (isnan(local_n[i])) {
	fprintf(stderr,"Population not a number local_n[%d]=%f\n",i,local_n[i]);
	errexit();
      }
      local_n[i] -= lbcoeff[i][0]*avg_rho;
      //fprintf(stderr,"local_n[%d]=%e\n",i,local_n[i]+lbcoeff[i][0]*avg_rho);      
    }
    else {
      local_n[i] = - lbcoeff[i][0]*avg_rho;
      //fprintf(stderr,"local_n[%d]=%e\n",i,local_n[i]+lbcoeff[i][0]*avg_rho);
    }
  }



  lb_calc_local_fields(node,1);
  LB_TRACE(fprintf(stderr,"Final: rho=%e\t",*(node->rho)));
  LB_TRACE(fprintf(stderr,"j=(%e,%f,%e)\t",node->j[0],node->j[1],node->j[2]));
  LB_TRACE(fprintf(stderr,"pi=(%e,%e,%e,%e,%e,%e)\n",node->pi[0],node->pi[1],node->pi[2],node->pi[3],node->pi[4],node->pi[5]));

  //if (fabs(node->pi[3])>ROUND_ERROR_PREC) {
  //  fprintf(stderr,"x-z-shear %e!\n",node->pi[3]);
  //  errexit();
  //}

  /* free memory */
  free(*A);
  free(A);
  free(b);
  free(perms);

}
#endif

MDINLINE void lb_set_boundary_node(int index, double rho, double *v, double *pi) {
#if 0 //problems with slip_pref (georg, 03.08.10)
  switch (lb_boundary_par.type) {

  case LB_BOUNDARY_NONE:
  case LB_BOUNDARY_BOUNCE_BACK:
  case LB_BOUNDARY_SPECULAR_REFLECTION:
  case LB_BOUNDARY_SLIP_REFLECTION:
    v[0] = v[1] = v[2] = 0.0;
    pi[0] = pi[1] = pi[2] = pi[3] = pi[4] = pi[5] = 0.0;
    lb_calc_n_equilibrium(index, 0.0, v, pi);
    break;
  case LB_BOUNDARY_PARTIAL_SLIP:
    lb_boundary_equilibrium(index, rho, v, pi, -1);
    break;

  }
#endif //if 0
}

/** Apply boundary conditions to the LB fluid. */
MDINLINE void lb_boundary_conditions() {
   lb_bounce_back();
#if 0 //problems with slip_pref (georg, 03.08.10)
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
    //lb_partial_slip();
    break;

  }
#endif //if 0
}


/** Parser for the \ref lbfluid command. */
int lbboundaries_cmd(ClientData data, Tcl_Interp *interp, int argc, char **argv);


extern int n_lb_boundaries;
extern LB_Boundary *lb_boundaries;

#endif /* LB_BOUNDARIES */

#endif /* LB_BOUNDARIES_H */
