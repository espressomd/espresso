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

/** \file lattice.h 
 *
 * Lattice data structures 
 * 
 */

#ifndef LATTICE_H
#define LATTICE_H

#include "utils.h"
#include "grid.h"

#ifdef LATTICE

/** Lattice off */
#define LATTICE_OFF  0

/** Lattice Boltzmann */
#define LATTICE_LB   1

/** Switch determining the type of lattice dynamics. A value of zero
 *  means that there is no lattice dynamics. Different types can be
 *  combined by or'ing the respective flags.
 *  So far, only \ref LATTICE_OFF and \ref LATTICE_LB exist.
 */
extern int lattice_switch ;

/** Data structure describing a lattice.
 *  Contains the lattice layout and pointers to the data fields.
 *  For parallelization purposes, it is assumed that a halo region
 *  surrounds the local lattice sites.
 */
typedef struct _Lattice {

  /** number of local lattice sites in each direction (excluding halo) */
  int grid[3] ;
  /** number of lattice sites in each direction including halo */
  int halo_grid[3] ;

  /** total number (volume) of local lattice sites (excluding halo) */
  int grid_volume ;
  /** total number (volume) of lattice sites (including halo) */
  int halo_grid_volume ;
  /** number of lattice sites in the halo region */
  int halo_grid_surface ;
  /** offset for number of halo sites stored in front of the local
      lattice sites */
  int halo_offset ;

  /** lattice constant */
  double agrid ;
  /** time step of lattice dynamics */
  double tau ;

  /** pointer to the fields living on the lattice.
   *  For complex lattice data, this can be a struct holding
   *  pointers to the actual data. */
  void *fields;
  /** pointer to the actual lattice data.
   *  This can be a contiguous field of arbitrary data. */
  void *data;

} Lattice;

/** Initialize lattice.
 *
 * This function initializes the variables describing the lattice
 * layout. Important: The lattice data is <em>not</em> allocated here!
 *
 * \param lattice pointer to the lattice
 * \param agrid   lattice spacing
 * \param tau     time step for lattice dynamics
 */
void init_lattice(Lattice *lattice, double agrid, double tau);

/** Map a global lattice site to the node grid.
 *
 *  This function determines the processor responsible for
 *  the specified lattice site. The coordinates of the site are
 *  taken as global coordinates and are returned as local coordinates.
 *
 * \param  lattice pointer to the lattice
 * \param  ind     global coordinates of the lattice site (Input)
 * \param  ind     local coordinates of the lattice site (Output)
 * \return         index of the node for the lattice site
 */
MDINLINE int map_lattice_to_node(Lattice *lattice, int *ind) {
  int im[3];
  
  /* determine coordinates in node_grid */
  im[0] = (int)floor(ind[0]*lattice->agrid*box_l_i[0]*node_grid[0]);
  im[1] = (int)floor(ind[1]*lattice->agrid*box_l_i[1]*node_grid[1]);
  im[2] = (int)floor(ind[2]*lattice->agrid*box_l_i[2]*node_grid[2]);

  /* change from global to local lattice coordinates */
  ind[0] -= im[0]*lattice->grid[0];
  ind[1] -= im[1]*lattice->grid[1];
  ind[2] -= im[2]*lattice->grid[2];

  /* return linear index into node array */
  return map_array_node(im);
}

/** Map a spatial position to the surrounding lattice sites.
 *
 * This function takes a global spatial position and determines the
 * surrounding elementary cell of the lattice for this position. 
 * The distance fraction in each direction is also calculated.
 * <br><em>Remarks:</em>
 * <ul>
 * <li>The spatial position has to be in the local domain.</li>
 * <li>The lattice sites of the elementary cell are returned as local indices</li>
 * </ul>
 * \param lattice    pointer to the lattice (Input)
 * \param pos        spatial position (Input)
 * \param node_index local indices of the surrounding lattice sites (Output)
 * \param delta      distance fraction of pos from the surrounding
 *                   elementary cell, 6 directions (Output)
 */
MDINLINE void map_position_to_lattice(Lattice *lattice, const double pos[3], int node_index[8], double delta[6]) {

  int x,y,z,dir,ind[3] ;
  double lpos, rel;

  /* determine the elementary lattice cell containing the particle
     and the relative position of the particle in this cell */ 
  for (dir=0;dir<3;dir++) {

    rel = (lpos = pos[dir] - my_left[dir])/lattice->agrid + 1.0; // +1 for halo offset
    ind[dir] = (int)floor(rel);
    
    /* surrounding elementary cell is not completely inside this box,
       adjust if this is due to round off errors */
    if (ind[dir] < 0) {
      fprintf(stderr,"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette.\n",this_node,pos[0],pos[1],pos[2]);
    }
    else if (ind[dir] > lattice->grid[dir]) {
      if (lpos - local_box_l[dir] < ROUND_ERROR_PREC)
	ind[dir] = lattice->grid[dir];
      else
	fprintf(stderr,"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette.\n",this_node,pos[0],pos[1],pos[2]);
    }

    delta[3+dir] = rel - ind[dir]; // delta_x/a
    delta[dir]   = 1.0 - delta[3+dir];   

  }

  for (x=0;x<2;x++) {
    for (y=0;y<2;y++) {
      for (z=0;z<2;z++) {
	node_index[z*4+y*2+x] = get_linear_index(ind[0]+x,ind[1]+y,ind[2]+z,lattice->halo_grid) ;
      }
    }
  }

}

#endif /* LATTICE */

#endif /* LATTICE_H */
