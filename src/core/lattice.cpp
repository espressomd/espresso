/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file lattice.cpp 
 *
 * Lattice class definition
 *
 */

#include "debug.hpp"
#include "lattice_inline.hpp"

int lattice_switch = LATTICE_OFF ;

int Lattice::init(double *agrid, double* offset, int halo_size, size_t dim) {
    /* determine the number of local lattice nodes */
    for (int d=0; d<3; d++) {
        this->agrid[d] = agrid[d];
        this->global_grid[d] = (int)dround(box_l[d]/agrid[d]);
        this->offset[d]=offset[d];
        this->local_index_offset[d]=(int) ceil((my_left[d]-this->offset[d])/this->agrid[d]);
        this->local_offset[d] = this->offset[d] +
            this->local_index_offset[d]*this->agrid[d];
        this->grid[d] = (int) ceil ( ( my_right[d] - this->local_offset[d]-ROUND_ERROR_PREC )
                                     / this->agrid[d]);
    }

    // sanity checks
    for (int dir=0;dir<3;dir++) {
      // check if local_box_l is compatible with lattice spacing
      if (fabs(local_box_l[dir]-this->grid[dir]*agrid[dir]) > ROUND_ERROR_PREC*box_l[dir]) {
          runtimeErrorMsg() << "Lattice spacing agrid["<< dir << "]=" << agrid[dir] \
              << " is incompatible with local_box_l["<< dir << "]=" << local_box_l[dir]\
              << " ( box_l["<< dir << "]=" << box_l[dir] \
              << " node_grid["<< dir << "]=" << node_grid[dir] <<" )";
      }
    }

    LATTICE_TRACE(fprintf(stderr,"%d: box_l (%.3f,%.3f,%.3f) grid (%d,%d,%d) node_neighbors (%d,%d,%d,%d,%d,%d)\n",this_node,local_box_l[0],local_box_l[1],local_box_l[2],this->grid[0],this->grid[1],this->grid[2],node_neighbors[0],node_neighbors[1],node_neighbors[2],node_neighbors[3],node_neighbors[4],node_neighbors[5]));

    this->halo_size = halo_size;
    /* determine the number of total nodes including halo */
    this->halo_grid[0] = this->grid[0] + 2*halo_size ;
    this->halo_grid[1] = this->grid[1] + 2*halo_size ;
    this->halo_grid[2] = this->grid[2] + 2*halo_size ;

    this->halo_grid_volume = this->halo_grid[0]*this->halo_grid[1]*this->halo_grid[2] ;
    this->halo_offset = get_linear_index(halo_size,halo_size,halo_size,this->halo_grid) ;

    return ES_OK;

}

void Lattice::map_position_to_lattice(const Vector3d& pos, index_t node_index[8], double delta[6]) {
    int ind[3];

    /* determine the elementary lattice cell containing the particle
       and the relative position of the particle in this cell */
    for (int dir=0;dir<3;dir++) {
        double lpos = pos[dir] - my_left[dir];
        double rel = lpos/this->agrid[dir] + 0.5; // +1 for halo offset
        ind[dir] = (int)floor(rel);

        /* surrounding elementary cell is not completely inside this box,
           adjust if this is due to round off errors */
        if (ind[dir] < 0) {
            if (fabs(rel) < ROUND_ERROR_PREC) {
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

/********************** static Functions **********************/

void Lattice::map_position_to_lattice_global (Vector3d& pos, int ind[3], double delta[6], double tmp_agrid) {
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
