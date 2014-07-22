/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
/** \file lattice.hpp 
 *
 * Lattice class definition
 * Contains the lattice layout and pointers to the data fields.
 * For parallelization purposes, it is assumed that a halo region
 * surrounds the local lattice sites.
 */

#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include <mpi.h>
#include "utils.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

/** Switch determining the type of lattice dynamics. A value of zero
 *  means that there is no lattice dynamics. Different types can be
 *  combined by or'ing the respective flags.
 *  So far, only \ref LATTICE_OFF and \ref LATTICE_LB exist.
 */
extern int lattice_switch;

#define LATTICE_LB   1  /** Lattice Boltzmann */
#define LATTICE_LB_GPU   2  /** Lattice Boltzmann */
#define INTERPOLATION_LINEAR 1
#define index_t long
#define LATTICE_OFF  0  /** Lattice off */

enum { LATTICE_ANISOTROPIC = 1, 
       LATTICE_X_NOTEXT = 2, LATTICE_Y_NOTEXT = 4, LATTICE_Z_NOTEXT = 8 };

class Lattice {
public:
  int grid[3] ;/** number of local lattice sites in each direction (excluding halo) */
  int global_grid[3];
  unsigned int dim;
  double agrid[3];/** lattice constant */

  int halo_grid[3] ;/** number of lattice sites in each direction including halo */
  int halo_size;/** halo size in all directions */

  double offset[3];/** global offset */
  double local_offset[3];
  int local_index_offset[3];

  unsigned int interpolation_type;
  char flags;
  size_t element_size;/** Size of each element in size units (=bytes) */
  size_t lattice_dim;/** Dimension of the field, assuming entries are arrays */

  index_t grid_volume;/** total number (volume) of local lattice sites (excluding halo) */
  index_t halo_grid_volume;/** total number (volume) of lattice sites (including halo) */
  index_t halo_grid_surface;/** number of lattice sites in the halo region */
  index_t halo_offset;/** offset for number of halo sites stored in front of the local lattice sites */

  void *_data;/** pointer to the actual lattice data. This can be a contiguous field of arbitrary data. */

  /** particle representation of this lattice. This is needed to
   *  specify interactions between particles and the lattice.
   *  Actually used are only the identity and the type. */
  Particle part_rep;

  /* Constructor */
  Lattice() {}

  /** Initialize lattice.
   *
   * This function initializes the variables describing the lattice
   * layout. Important: The lattice data is <em>not</em> allocated here!
   *
   * \param lattice pointer to the lattice
   * \param agrid   lattice spacing
   */
  int init(double* agrid, double* offset, int halo_size, size_t element_size);

  /** lattice memory allocation.
   * \param lattice pointer to the lattice
   */
  void allocate_memory();
  void allocate_memory(size_t element_size);

  void interpolate(double* pos, double* value);

  void interpolate_linear(double* pos, double* value);

  void interpolate_gradient(double* pos, double* value);

  void interpolate_linear_gradient(double* pos, double* value);

  void set_data_for_global_position_with_periodic_image(double* pos, void* data);

  /********************** Inline Functions **********************/

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
  int map_lattice_to_node(int *ind, int *grid) {

    /* determine coordinates in node_grid */
    grid[0] = (int)floor(ind[0]*this->agrid[0]*box_l_i[0]*node_grid[0]);
    grid[1] = (int)floor(ind[1]*this->agrid[1]*box_l_i[1]*node_grid[1]);
    grid[2] = (int)floor(ind[2]*this->agrid[2]*box_l_i[2]*node_grid[2]);

    //fprintf(stderr,"%d: (%d,%d,%d)\n",this_node,grid[0],grid[1],grid[2]);

    /* change from global to local lattice coordinates */
    ind[0] = ind[0] - grid[0]*this->grid[0] + this->halo_size;
    ind[1] = ind[1] - grid[1]*this->grid[1] + this->halo_size;
    ind[2] = ind[2] - grid[2]*this->grid[2] + this->halo_size;

    /* return linear index into node array */
    return map_array_node(grid);
  }

  int global_pos_in_local_box(double pos[3]) {
    if (!(pos[0]>my_left[0]  &&  pos[0]<my_right[0] ))
      return 0;
    if (!(pos[1]>my_left[1]  &&  pos[1]<my_right[1] ))
      return 0;
    if (!(pos[2]>my_left[2]  &&  pos[2]<my_right[2] ))
      return 0;
    return 1;
  }

  int global_pos_in_local_halobox(double pos[3]) {
    for (unsigned int i=0; i<this->dim; i++) {
      if (!(pos[i]>this->local_offset[i]-this->halo_size*this->agrid[i] &&
            pos[i]<this->local_offset[i]+this->halo_grid[i]*this->agrid[i] ))
        return 0;
    }
    return 1;
  }

  int global_pos_to_lattice_index_checked(double pos[3], int* index) {
    int i;
    for (i=0; i<3; i++)
      if (abs(fmod(pos[i]-this->offset[i],this->agrid[i])) > ROUND_ERROR_PREC)
        return ES_ERROR;
    int ind[3];
    for (i=0; i<3; i++)
      ind[i] = (int) round((pos[i]-this->offset[i])/this->agrid[i]);
    *index = get_linear_index(this->halo_size + ind[0], this->halo_size + ind[1], this->halo_size + ind[2], this->halo_grid);
    return ES_OK;
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
  /* @TODO: Implement! */
  int map_lattice_to_position(int *ind, int *grid) {
    return 0;
  }

  void map_linear_index_to_global_pos(index_t ind, double* pos) {
    int index_in_halogrid[3];
    get_grid_pos(ind, &index_in_halogrid[0], &index_in_halogrid[1], &index_in_halogrid[2], this->halo_grid);
    pos[0] = this->local_offset[0] + (index_in_halogrid[0] - this->halo_size)*this->agrid[0];
    pos[1] = this->local_offset[1] + (index_in_halogrid[1] - this->halo_size)*this->agrid[1];
    pos[2] = this->local_offset[2] + (index_in_halogrid[2] - this->halo_size)*this->agrid[2];
  }

  void map_local_index_to_pos(index_t* index, double* pos) {
    pos[0] = this->local_offset[0] + (index[0])*this->agrid[0];
    pos[1] = this->local_offset[1] + (index[1])*this->agrid[1];
    pos[2] = this->local_offset[2] + (index[2])*this->agrid[2];
  }

  int map_global_index_to_halo_index(index_t* global_index, index_t* halo_index) {
    int out=0;
    for (int d=0; d<3; d++) {
      halo_index[d] = global_index[d]-this->local_index_offset[d] +this->halo_size;
      if (halo_index[d] < 0 || halo_index[d] >= this->halo_grid[d])
        out=1;
    }

    if (out) {
      return 0;
    }
    return 1;

  }

  void map_halo_index_to_pos(index_t* index_in_halogrid, double* pos) {
    pos[0] = this->local_offset[0] + (index_in_halogrid[0] - this->halo_size)*this->agrid[0];
    pos[1] = this->local_offset[1] + (index_in_halogrid[1] - this->halo_size)*this->agrid[1];
    pos[2] = this->local_offset[2] + (index_in_halogrid[2] - this->halo_size)*this->agrid[2];
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
  void map_position_to_lattice(const double pos[3], index_t node_index[8], double delta[6]) {

    int dir,ind[3] ;
    double lpos, rel;

    /* determine the elementary lattice cell containing the particle
       and the relative position of the particle in this cell */
    for (dir=0;dir<3;dir++) {
      lpos = pos[dir] - my_left[dir];
      rel = lpos/this->agrid[dir] + 0.5; // +1 for halo offset
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
      else if (ind[dir] > this->grid[dir]) {
        if (lpos - local_box_l[dir] < ROUND_ERROR_PREC*local_box_l[dir])
      ind[dir] = this->grid[dir];
        else
      fprintf(stderr,"%d: map_position_to_lattice: position (%f,%f,%f) not inside a local plaquette in dir %d ind[dir]=%d rel=%f lpos=%f.\n",this_node,pos[0],pos[1],pos[2],dir,ind[dir],rel,lpos);
      }

      delta[3+dir] = rel - ind[dir]; // delta_x/a
      delta[dir]   = 1.0 - delta[3+dir];
    }

    node_index[0] = get_linear_index(ind[0],ind[1],ind[2],this->halo_grid);
    node_index[1] = node_index[0] + 1;
    node_index[2] = node_index[0] + this->halo_grid[0];
    node_index[3] = node_index[0] + this->halo_grid[0] + 1;
    node_index[4] = node_index[0] + this->halo_grid[0]*this->halo_grid[1];
    node_index[5] = node_index[4] + 1;
    node_index[6] = node_index[4] + this->halo_grid[0];
    node_index[7] = node_index[4] + this->halo_grid[0] + 1;
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
  static void map_position_to_lattice_global (double pos[3], int ind[3], double delta[6], double tmp_agrid) {
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


  void get_data_for_halo_index(index_t* ind, void** data) {
    (*data) = ((char*)this->_data) + get_linear_index(ind[0], ind[1], ind[2], this->halo_grid)*this->element_size;
  }

  void get_data_for_linear_index(index_t ind, void** data) {
    (*data) = ((char*)this->_data) + ind*this->element_size;
  }

  void get_data_for_local_index(index_t* ind, void** data) {
    index_t index_in_halogrid[3];
    index_in_halogrid[0] = ind[0]+this->halo_size;
    index_in_halogrid[1] = ind[1]+this->halo_size;
    index_in_halogrid[2] = ind[2]+this->halo_size;
    (*data) = ((char*)this->_data) + get_linear_index(index_in_halogrid[0], index_in_halogrid[1], index_in_halogrid[2], this->halo_grid)*this->element_size;
  }

  void set_data_for_local_halo_grid_index(index_t* ind, void* data) {
    memcpy(((char*)this->_data) + get_linear_index(ind[0], ind[1], ind[2],  this->halo_grid)*this->element_size, data, this->element_size);

  }

  void set_data_for_local_grid_index(index_t* ind, void* data) {
    memcpy(((char*)this->_data) + get_linear_index(ind[0]+this->halo_size, ind[1]+this->halo_size, ind[2]+this->halo_size,  this->halo_grid)*this->element_size, data, this->element_size);
  }

  int global_pos_to_lattice_halo_index(double* pos, index_t*  ind) {
    for (int i = 0; i<3; i++) {
      ind[i] = (int) floor((pos[i]-this->local_offset[i])/this->agrid[i]+ROUND_ERROR_PREC)+this->halo_size;
      if (ind[i] < 0 || ind[i] >= this->halo_grid[i])
        return 0;
    }
    return 1;
  }

};

#endif /* LATTICE_HPP */
