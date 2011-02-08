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

/** \file lb_boundaries_gpu.h
 *
 * Boundary conditions for Lattice Boltzmann GPU fluid dynamics.
 * Header file for \ref lb_boundaries_gpu.c.
 *
 */

#ifndef LB_BOUNDARIES_GPU_H
#define LB_BOUNDARIES_GPU_H
#include <tcl.h>
#include "utils.h"
//#include "halo.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "constraint.h"
#ifdef __cplusplus
}
#endif
#include "config.h"
#include "lbgpu.h"

#ifdef LB_BOUNDARIES_GPU

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
} LB_boundary_gpu;

/*@}*/

#ifdef __cplusplus
extern "C" {
#endif
void lb_init_boundaries_gpu();

//void lb_init_boundaries_GPU(int host_number_of_boundnodes, int *host_boundindex);

#ifdef __cplusplus
}
#endif

/*@}*/

int tclcommand_lbboundary_gpu(Tcl_Interp *interp, int argc, char **argv);
#endif /* LB_BOUNDARIES_GPU */
#endif /* LB_BOUNDARIES_GPU_H */

