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

LB_Boundary lb_boundary_par = { LB_BOUNDARY_NONE, 0.0 };

/** Initialize a planar boundary specified by a wall constraint.
 * @param plane The \ref Constraint_wall struct describing the boundary.
 */
static void lb_init_constraint_wall(Constraint_wall* plane) {

  int x, y, z, i, index, mask;
  double pos[3], dist;

  mask = (plane->d < 0.0) ? -1 : +1;

  for (x=0;x<lblattice.halo_grid[0];x++) {
    for (y=0;y<lblattice.halo_grid[1];y++) {
      for (z=0;z<lblattice.halo_grid[2];z++) {

	pos[0] = my_left[0] + (x-1)*lblattice.agrid;
	pos[1] = my_left[1] + (y-1)*lblattice.agrid;
	pos[2] = my_left[2] + (z-1)*lblattice.agrid;

	dist = scalar(pos,plane->n) - plane->d;

	if (fabs(dist) < lblattice.agrid - ROUND_ERROR_PREC) {
	  printf("%d: (%2d, %2d, %2d) = %e (%d) (%f,%f,%f)\n",this_node,x,y,z,fabs(dist),mask,pos[0],pos[1],pos[2]);
	  index = get_linear_index(x,y,z,lblattice.halo_grid);
	  lbfluid[index].boundary = mask;
	  lbfluid[index].nvec     = plane->n;
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


static int lbboundaries_parse_slip_reflection(Tcl_Interp *interp, int argc, char **argv) {

  if (argc <1) {
    Tcl_AppendResult(interp, "lbboundaries slip_reflection requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(lb_boundary_par.slip_pref)) {
    Tcl_AppendResult(interp, "wrong argument for lbboundaries slip_reflection", (char *)NULL);
    return TCL_ERROR;
  }
  
  lb_boundary_par.type = LB_BOUNDARY_SLIP_REFLECTION;

  return TCL_OK;

}

static int lbboundaries_parse_partial_slip(Tcl_Interp *interp, int argc, char **argv) {

  if (argc <1) {
    Tcl_AppendResult(interp, "lbboundaries slip requires 1 argument", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(lb_boundary_par.slip_pref)) {
    Tcl_AppendResult(interp, "wrong argument for lbboundaries slip", (char *)NULL);
    return TCL_ERROR;
  }
  
  lb_boundary_par.type = LB_BOUNDARY_PARTIAL_SLIP;

  return TCL_OK;

}

/** Parser for the \ref lbfluid command. */
int lbboundaries_cmd(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int err = TCL_ERROR;

  if (argc < 2) {
    Tcl_AppendResult(interp, "too few arguments to \"lbboundaries\"", (char *)NULL);
    err = TCL_ERROR;
  }
  else if (ARG1_IS_S("off")) {
    err = TCL_ERROR;
  }
  else if (ARG1_IS_S("bounce_back")) {
    lb_boundary_par.type = LB_BOUNDARY_BOUNCE_BACK;
    err = TCL_OK;
  }
  else if (ARG1_IS_S("specular_reflections")) {
    lb_boundary_par.type = LB_BOUNDARY_SPECULAR_REFLECTION;
    err = TCL_OK;
  }
  else if (ARG1_IS_S("slip_reflection")) {
    err = lbboundaries_parse_slip_reflection(interp, argc-2, argv+2);
  }
  else if (ARG1_IS_S("partial_slip")) {
    err = lbboundaries_parse_partial_slip(interp, argc-2, argv+2);
  }
  else {
    Tcl_AppendResult(interp, "unkown boundary condition \"", argv[1], (char *)NULL);
    err = TCL_ERROR;
  }

  return err;

}

#endif /* CONSTRAINTS */
#endif /* LB */
