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

/** \file lb-boundaries.c
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.h.
 *
 */

#include "utils.h"
#include "constraint.h"
#include "lb-boundaries.h"

#ifdef LB
#ifdef CONSTRAINTS

/** Initialize a planar boundary specified by a wall constraint.
 * @param plane The \ref Constraint_wall struct describing the boundary.
 */
static void lb_init_constraint_wall(Constraint_wall* plane) {

  int x, y, z;
  double pos[3], dist;

  for (x=0;x<lblattice.halo_grid[0];x++) {
    for (y=0;y<lblattice.halo_grid[1];y++) {
      for (z=0;z<lblattice.halo_grid[2];z++) {

	pos[0] = my_left[0] + (x-1)*lblattice.agrid;
	pos[1] = my_left[1] + (y-1)*lblattice.agrid;
	pos[2] = my_left[2] + (z-1)*lblattice.agrid;

	dist = scalar(pos,plane->n) - plane->d;

	if (fabs(dist) < lblattice.agrid) {
	  //printf("%d %d %d\n",x,y,z);
	  lbfluid[get_linear_index(x,y,z,lblattice.halo_grid)].boundary = 1;
	}

      }
    }
  }

}

/** Initialize boundary conditions for all constraints in the system. */
void lb_init_constraints() {

  int n;
  char *errtxt;

  for (n=0;n<lblattice.halo_grid_volume;n++) {
    lbfluid[n].boundary = 0;
  }

  for (n=0;n<n_constraints;n++) {
    switch (constraints[n].type) {
    case CONSTRAINT_WAL:
      lb_init_constraint_wall(&constraints[n].c.wal);
      break;
    default:
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{109 constraint type %d not implemented in lb_init_constraints()\n",constraints[n].type);
    }
  }

}

#endif /* CONSTRAINTS */
#endif /* LB */
