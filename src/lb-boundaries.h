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
/** \file lb-boundaries.h
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.c.
 *
 * In the current version only simple bounce back walls are implemented. Thus
 * after the streaming step, in all wall nodes all populations are bounced
 * back from where they came from. Ulf Schiller spent a lot of time
 * working on more powerful alternatives, they are to be found in the
 * lb_testing branch of espresso until the end of 2010. Now we stripped
 * down the code to a minimum, as most of it was not sufficiently understandable.
 *
 * Anyone who wants to revive these, please look into the git.
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

/** wall constraint applied */
#define LB_BOUNDARY_WAL 1
/** spherical constraint applied */
#define LB_BOUNDARY_SPH 2
/** (finite) cylinder shaped constraint applied */
#define LB_BOUNDARY_CYL 3
/** a pore geometry */
#define LB_BOUNDARY_POR 4

// If we have several possible types of boundary treatment
#define LB_BOUNDARY_BOUNCE_BACK 1

/** Parser for the \ref lbfluid command. */
int tclcommand_lbboundary(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Structure to specify a boundary. */
typedef struct {
  /** type of the boundary. */
  int type;
  double slip_pref;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
    Constraint_pore pore;
  } c;
  double force[3];
  double velocity[3];
} LB_Boundary;

extern int n_lb_boundaries;
extern LB_Boundary *lb_boundaries;

/*@}*/

/** Initializes the constrains in the system. 
 *  This function determines the lattice sited which belong to boundaries
 *  and marks them with a corresponding flag. 
 */
void lb_init_boundaries();
#endif // LB_BOUNDARIES
int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv);
#ifdef LB_BOUNDARIES

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
  int k,i,l;
  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0]*lblattice.halo_grid[1];
  int next[19];
  double population_shift;
  double rho, j[3], pi[6];
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
      lb_calc_local_fields(k, &rho, j, pi);
      lb_boundaries[lbfields[k].boundary-1].force[0]+=2*j[0];
      lb_boundaries[lbfields[k].boundary-1].force[1]+=2*j[1];
      lb_boundaries[lbfields[k].boundary-1].force[2]+=2*j[2];
      for (i=0; i<19; i++) {
        population_shift=0;
        for (l=0; l<3; l++)
          population_shift-=lbpar.agrid*lbpar.agrid*lbpar.agrid*lbpar.rho*2*lbmodel.c[i][l]*lb_boundaries[lbfields[k].boundary-1].velocity[l];
        lbfluid[1][reverse[i]][k-next[i]]   = lbfluid[1][i][k] + population_shift;
      }
    }
  }
  
#else
#error Bounce back boundary conditions are only implemented for PUSH scheme!
#endif
#else
#error Bounce back boundary conditions are only implemented for D3Q19!
#endif
}



#endif /* LB_BOUNDARIES */

#endif /* LB_BOUNDARIES_H */
