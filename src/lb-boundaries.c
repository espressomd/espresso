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
