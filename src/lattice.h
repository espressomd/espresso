/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file lattice.h 
 *
 * Lattice data structures 
 * 
 */
/** Lattice Boltzmann */
#define LATTICE_LB   1

/** Lattice Boltzmann */
#define LATTICE_LB_GPU   2

#ifndef LATTICE_H
#define LATTICE_H

#include <mpi.h>
#include "utils.h"
#include "grid.h"
#include "particle_data.h"

extern int lattice_switch;

#define index_t long

/** Lattice off */
#define LATTICE_OFF  0

#ifdef LATTICE

/** Lattice Boltzmann */
//#define LATTICE_LB   1

/** Lattice Boltzmann */
//#define LATTICE_LB_GPU   2

/** Switch determining the type of lattice dynamics. A value of zero
 *  means that there is no lattice dynamics. Different types can be
 *  combined by or'ing the respective flags.
 *  So far, only \ref LATTICE_OFF and \ref LATTICE_LB exist.
 */

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
  index_t grid_volume;
  /** total number (volume) of lattice sites (including halo) */
  index_t halo_grid_volume;
  /** number of lattice sites in the halo region */
  index_t halo_grid_surface;
  /** offset for number of halo sites stored in front of the local
      lattice sites */
  index_t halo_offset;

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

  /** particle representation of this lattice. This is needed to
   *  specify interactions between particles and the lattice.
   *  Actually used are only the identity and the type. */
  Particle part_rep;

  /** datatypes that describe the data layout of the lattice. */
  MPI_Datatype datatype;
  struct _Fieldtype *fieldtype;

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
 * \param  grid     local coordinates of the lattice site (Output)
 * \return         index of the node for the lattice site
 */
MDINLINE int map_lattice_to_node(Lattice *lattice, int *ind, int *grid) {
  
  /* determine coordinates in node_grid */
  grid[0] = (int)floor(ind[0]*lattice->agrid*box_l_i[0]*node_grid[0]);
  grid[1] = (int)floor(ind[1]*lattice->agrid*box_l_i[1]*node_grid[1]);
  grid[2] = (int)floor(ind[2]*lattice->agrid*box_l_i[2]*node_grid[2]);

  //fprintf(stderr,"%d: (%d,%d,%d)\n",this_node,grid[0],grid[1],grid[2]);

  /* change from global to local lattice coordinates (+1 is for halo) */
  ind[0] = ind[0] - grid[0]*lattice->grid[0] + 1;
  ind[1] = ind[1] - grid[1]*lattice->grid[1] + 1;
  ind[2] = ind[2] - grid[2]*lattice->grid[2] + 1;

  /* return linear index into node array */
  return map_array_node(grid);
}

/** Map a local lattice site to the global position.
 *
 *  This function determines the processor responsible for
 *  the specified lattice site. The coordinates of the site are
 *  taken as global coordinates andcoordinates of the site are
 *  taken as global coordinates and are returned as local coordinates.
 *
 * \param  lattice pointer to the lattice
 * \param  ind     global coordinates of the lattice site (Input)
 * \param  grid    local coordinates of the lattice site (Output)
 * \return         index of the node for the lattice site
 */
MDINLINE int map_lattice_to_position(Lattice *lattice, int *ind, int *grid) {
  return 0;
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
MDINLINE void map_position_to_lattice(Lattice *lattice, const double pos[3], index_t node_index[8], double delta[6]) {

  int dir,ind[3] ;
  double lpos, rel;

  /* determine the elementary lattice cell containing the particle
     and the relative position of the particle in this cell */ 
  for (dir=0;dir<3;dir++) {

    rel = (lpos = pos[dir] - my_left[dir])/lattice->agrid + 0.5; // +1 for halo offset
    ind[dir] = (int)floor(rel);
    
    /* surrounding elementary cell is not completely inside this box,
       adjust if this is due to round off errors */
    if (ind[dir] < 0) {
      if (abs(rel) < ROUND_ERROR_PREC) {
	ind[dir] = 0; // TODO
      } else {
	fprintf(stderr,"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette in dir %d ind[dir]=%d rel=%f lpos=%f.\n",this_node,pos[0],pos[1],pos[2],dir,ind[dir],rel,lpos);
      }
    }
    else if (ind[dir] > lattice->grid[dir]) {
      if (lpos - local_box_l[dir] < ROUND_ERROR_PREC*local_box_l[dir])
	ind[dir] = lattice->grid[dir];
      else
	fprintf(stderr,"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette in dir %d ind[dir]=%d rel=%f lpos=%f.\n",this_node,pos[0],pos[1],pos[2],dir,ind[dir],rel,lpos);
    }

    delta[3+dir] = rel - ind[dir]; // delta_x/a
    delta[dir]   = 1.0 - delta[3+dir];   

  }

  node_index[0] = get_linear_index(ind[0],ind[1],ind[2],lattice->halo_grid);
  node_index[1] = node_index[0] + 1;
  node_index[2] = node_index[0] + lattice->halo_grid[0];
  node_index[3] = node_index[0] + lattice->halo_grid[0] + 1;
  node_index[4] = node_index[0] + lattice->halo_grid[0]*lattice->halo_grid[1];
  node_index[5] = node_index[4] + 1;
  node_index[6] = node_index[4] + lattice->halo_grid[0];
  node_index[7] = node_index[4] + lattice->halo_grid[0] + 1;

}

/** Map a spatial position to the surrounding lattice sites.
 *
 * This function takes a global spatial position and determines the
 * surrounding elementary cell of the lattice for this position. 
 * The distance fraction in each direction is also calculated.
 * <br><em>Remarks:</em>
 * <ul>
 * <li>The spatial position is given in global coordinates.</li>
 * <li>The lattice sites of the elementary cell are returned as local indices</li>
 * </ul>
 * \param pos        spatial position (Input)
 * \param ind        global index of the lower left lattice site (Output)
 * \param delta      distance fraction of pos from the surrounding
 *                   elementary cell, 6 directions (Output)
 * \param tmp_agrid  lattice mesh distance
 */

MDINLINE void map_position_to_lattice_global (double pos[3], int ind[3], double delta[6], double tmp_agrid) {
//not sure why I don't have access to agrid here so I make a temp var and pass it to this function
  int i;
  double rel[3];
  // fold the position onto the local box, note here ind is used as a dummy variable
	for (i=0;i<3;i++) {
		pos[i] = pos[i]-0.5*tmp_agrid;
	}
	
  fold_position (pos,ind);

  // convert the position into lower left grid point
	for (i=0;i<3;i++) {
		rel[i] = (pos[i])/tmp_agrid;
	}

  // calculate the index of the position
  for (i=0;i<3;i++) {
    ind[i] = floor(rel[i]);
  }

  // calculate the linear interpolation weighting
  for (i=0;i<3;i++) {
    delta[3+i] = rel[i] - ind[i];
    delta[i] = 1 - delta[3+i];
  }

}

#endif /* LATTICE */

#endif /* LATTICE_H */
