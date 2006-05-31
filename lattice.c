/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
 * and by which you are legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * You should have received a copy of that license along with this program;
 * if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
 * write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
 * Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
 */

/** \file lattice.c 
 *
 * Lattice data structures
 *
 */

#include "utils.h"
#include "grid.h"
#include "lattice.h"

#ifdef LATTICE

/** Switch determining the type of lattice dynamics. A value of zero
 *  means that there is no lattice dynamics. Different types can be
 *  combined by or'ing the respective flags.
 *  So far, only \ref LATTICE_OFF and \ref LATTICE_LB exist.
 */
int lattice_switch = LATTICE_OFF ;

/** Initialize lattice.
 *
 * This function initializes the variables describing the lattice
 * layout. Important: The lattice data is <em>not</em> allocated here!
 *
 * \param lattice pointer to the lattice
 * \param agrid   lattice spacing
 * \param tau     time step for lattice dynamics
 */
void init_lattice(Lattice *lattice, double agrid, double tau) {

  int dir;

  /* determine the number of local lattice nodes */
  lattice->grid[0] = local_box_l[0]/agrid;
  lattice->grid[1] = local_box_l[1]/agrid;
  lattice->grid[2] = local_box_l[2]/agrid;

  /* sanity checks */
  for (dir=0;dir<3;dir++) {
    /* check if local_box_l is compatible with lattice spacing */
    if (fabs(local_box_l[dir]-lattice->grid[dir]*agrid) > ROUND_ERROR_PREC) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{097 Lattice spacing agrid=%f is incompatible with local_box_l[%d]=%f (box_l[%d]=%f node_grid[%d]=%d) %f} ",agrid,dir,local_box_l[dir],dir,box_l[dir],dir,node_grid[dir],local_box_l[dir]-lattice->grid[dir]*agrid);
      return;
    }
  }

  /* set the lattice spacing */
  lattice->agrid = agrid ;
  lattice->tau = tau ;

  LATTICE_TRACE(fprintf(stderr,"%d: box_l (%.3f,%.3f,%.3f) grid (%d,%d,%d) node_neighbors (%d,%d,%d,%d,%d,%d)\n",this_node,local_box_l[0],local_box_l[1],local_box_l[2],lattice->grid[0],lattice->grid[1],lattice->grid[2],node_neighbors[0],node_neighbors[1],node_neighbors[2],node_neighbors[3],node_neighbors[4],node_neighbors[5]));

  /* determine the number of total nodes including halo */
  lattice->halo_grid[0] = lattice->grid[0] + 2 ;
  lattice->halo_grid[1] = lattice->grid[1] + 2 ;
  lattice->halo_grid[2] = lattice->grid[2] + 2 ;

  lattice->grid_volume = lattice->grid[0]*lattice->grid[1]*lattice->grid[2] ;
  lattice->halo_grid_volume = lattice->halo_grid[0]*lattice->halo_grid[1]*lattice->halo_grid[2] ;
  lattice->halo_grid_surface = lattice->halo_grid_volume - lattice->grid_volume ;
  lattice->halo_offset = get_linear_index(1,1,1,lattice->halo_grid) ;

}

#endif /* LATTICE */
