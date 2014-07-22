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
/** \file lattice.cpp 
 *
 * Lattice class definition
 *
 */

#include "utils.hpp"
#include "grid.hpp"
#include "lattice.hpp"

int lattice_switch = LATTICE_OFF ;

int Lattice::init(double *agrid, double* offset, int halo_size, size_t dim) {
  int dir;

  this->dim=dim;

 /* determine the number of local lattice nodes */
  for (int d=0; d<3; d++) {
    this->agrid[d] = agrid[d];
    this->global_grid[d] = (int)floor(box_l[d]/agrid[d]+ROUND_ERROR_PREC);
    this->offset[d]=offset[d];
    this->local_index_offset[d]=(int) ceil((my_left[d]-this->offset[d])/this->agrid[d]);
    this->local_offset[d] = this->offset[d] +
      this->local_index_offset[d]*this->agrid[d];
    this->grid[d] = (int) ceil ( ( my_right[d] - this->local_offset[d]-ROUND_ERROR_PREC )
      / this->agrid[d]);
  }

  /* sanity checks */
  for (dir=0;dir<3;dir++) {
    /* check if local_box_l is compatible with lattice spacing */
//    if (fabs(local_box_l[dir]-this->grid[dir]*agrid[dir]) > ROUND_ERROR_PREC*box_l[dir]) {
//      char *errtxt = runtime_error(128);
//      ERROR_SPRINTF(errtxt, "{097 Lattice spacing agrid[%d]=%f is incompatible with local_box_l[%d]=%f (box_l[%d]=%f node_grid[%d]=%d) %f} ",dir,agrid[dir],dir,local_box_l[dir],dir,box_l[dir],dir,node_grid[dir],local_box_l[dir]-this->grid[dir]*agrid[dir]);
//      return ES_ERROR;
//    }
  /* set the lattice spacing */
  }

  this->element_size = this->dim*sizeof(double);

  LATTICE_TRACE(fprintf(stderr,"%d: box_l (%.3f,%.3f,%.3f) grid (%d,%d,%d) node_neighbors (%d,%d,%d,%d,%d,%d)\n",this_node,local_box_l[0],local_box_l[1],local_box_l[2],this->grid[0],this->grid[1],this->grid[2],node_neighbors[0],node_neighbors[1],node_neighbors[2],node_neighbors[3],node_neighbors[4],node_neighbors[5]));

  this->halo_size = halo_size;
  /* determine the number of total nodes including halo */
  this->halo_grid[0] = this->grid[0] + 2*halo_size ;
  this->halo_grid[1] = this->grid[1] + 2*halo_size ;
  this->halo_grid[2] = this->grid[2] + 2*halo_size ;

  this->grid_volume = this->grid[0]*this->grid[1]*this->grid[2] ;
  this->halo_grid_volume = this->halo_grid[0]*this->halo_grid[1]*this->halo_grid[2] ;
  this->halo_grid_surface = this->halo_grid_volume - this->grid_volume ;
  this->halo_offset = get_linear_index(halo_size,halo_size,halo_size,this->halo_grid) ;

  this->interpolation_type = INTERPOLATION_LINEAR;

  allocate_memory();
  return ES_OK;

}

void Lattice::allocate_memory() {

  this->_data = malloc(this->element_size*this->halo_grid_volume);
  memset(this->_data, (unsigned int)(-1), this->element_size*this->halo_grid_volume);

}

void Lattice::interpolate(double* pos, double* value) {
  if (this->interpolation_type == INTERPOLATION_LINEAR) {
    interpolate_linear(pos, value);
  } else {
    char* c = runtime_error(128);
    ERROR_SPRINTF(c, "Unknown interpolation type");
  }
}

void Lattice::interpolate_linear(double* pos, double* value) {
  int left_halo_index[3];
  double d[3];
   if (this->halo_size <= 0) {
     char* c = runtime_error(128);
     ERROR_SPRINTF(c, "Error in interpolate_linear: halo size is 0");
     return;
   }
   for (int dim = 0; dim<3; dim++) {
     left_halo_index[dim]=(int) floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]) + this->halo_size;
     d[dim]=((pos[dim]-this->local_offset[dim])/this->agrid[dim] - floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]));
     if (left_halo_index[dim] < 0 || left_halo_index[dim] >= this->halo_grid[dim]) {
       char* c = runtime_error(128);
       ERROR_SPRINTF(c, "Error in interpolate_linear: Particle out of range");
       return;
     }
   }
   double w[8];
   index_t index[8];
   w[0] = (1-d[0])*(1-d[1])*(1-d[2]);
   index[0]=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2], this->halo_grid);
   w[1] = ( +d[0])*(1-d[1])*(1-d[2]);
   index[1]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2], this->halo_grid);
   w[2] = (1-d[0])*( +d[1])*(1-d[2]);
   index[2]=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2], this->halo_grid);
   w[3] = ( +d[0])*( +d[1])*(1-d[2]);
   index[3]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2], this->halo_grid);

   w[4] = (1-d[0])*(1-d[1])*( +d[2]);
   index[4]=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2]+1, this->halo_grid);
   w[5] = ( +d[0])*(1-d[1])*( +d[2]);
   index[5]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2]+1, this->halo_grid);
   w[6] = (1-d[0])*( +d[1])*( +d[2]);
   index[6]=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);
   w[7] = ( +d[0])*( +d[1])*( +d[2]);
   index[7]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);

   for (unsigned int i = 0; i<this->dim; i++) {
     value[i] = 0;
   }

   double* local_value;
   for (unsigned int i=0; i<8; i++) {
     get_data_for_linear_index(index[i], (void**) &local_value);
     for (unsigned int j = 0; j<this->dim; j++) {
       value[j]+=w[i]*local_value[j];
     }
   }
}

void Lattice::interpolate_gradient(double* pos, double* value) {
  if (this->interpolation_type == INTERPOLATION_LINEAR) {
    interpolate_linear_gradient(pos, value);
  } else {
    char* c = runtime_error(128);
    ERROR_SPRINTF(c, "Unknown interpolation type");
  }
}

void Lattice::interpolate_linear_gradient(double* pos, double* value) {
  int left_halo_index[3];
  double d[3];
   if (this->halo_size <= 0) {
     char* c = runtime_error(128);
     ERROR_SPRINTF(c, "Error in interpolate_linear: halo size is 0");
     return;
   }
   for (int dim = 0; dim<3; dim++) {
     left_halo_index[dim]=(int) floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]) + this->halo_size;
     d[dim]=((pos[dim]-this->local_offset[dim])/this->agrid[dim] - floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]));
     if (left_halo_index[dim] < 0 || left_halo_index[dim] >= this->halo_grid[dim]) {
       char* c = runtime_error(128);
       ERROR_SPRINTF(c, "Error in interpolate_linear: Particle out of range");
       return;
     }
   }
   
   index_t index;
   double* local_value;

   for (unsigned int i = 0; i<3*this->dim; i++) {
     value[i] = 0;
   }
   
   index=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2], this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*(1-d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= (1-d[0])*( -1   )*(1-d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= (1-d[0])*(1-d[1])*(  -1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2], this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*(1-d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= ( +d[0])*( -1   )*(1-d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= ( +d[0])*(1-d[1])*(  -1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2], this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*( +d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= (1-d[0])*( +1   )*(1-d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= (1-d[0])*( +d[1])*(  -1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2], this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*( +d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= ( +d[0])*( +1   )*(1-d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= ( +d[0])*( +d[1])*(  -1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]  , left_halo_index[1]  , left_halo_index[2] + 1, this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*(1-d[1])*( +d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= (1-d[0])*( -1   )*( +d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= (1-d[0])*(1-d[1])*(  +1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2]+1, this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*(1-d[1])*( +d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= ( +d[0])*( -1   )*( +d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= ( +d[0])*(1-d[1])*(  +1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*( +d[1])*( +d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= (1-d[0])*( +1   )*( +d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= (1-d[0])*( +d[1])*(  +1  ) * local_value[i] / this->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);
   for (unsigned int i = 0; i<this->dim; i++) {
     get_data_for_linear_index(index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*( +d[1])*( +d[2]) * local_value[i] / this->agrid[0];
     value[3*i+1]+= ( +d[0])*( +1   )*( +d[2]) * local_value[i] / this->agrid[1];
     value[3*i+2]+= ( +d[0])*( +d[1])*(  +1  ) * local_value[i] / this->agrid[2];
   }

}

void Lattice::set_data_for_global_position_with_periodic_image(double* pos, void* data) {

  index_t replica[3];
  index_t global_index[3];
  index_t halo_index[3];


  for (int i = 0; i<3; i++) {
    global_index[i] = (int) floor((pos[i]-this->offset[i])/this->agrid[i]+ROUND_ERROR_PREC);
  }

  for (int i=-2; i<=2; i++) {
    for (int j=-2; j<=2; j++) {
      for (int k=-2; k<=2; k++) {
        replica[0]=global_index[0]+i*this->global_grid[0];
        replica[1]=global_index[1]+j*this->global_grid[1];
        replica[2]=global_index[2]+k*this->global_grid[2];
        if (map_global_index_to_halo_index(replica, halo_index)) {
          set_data_for_local_halo_grid_index(halo_index, data);
        }

      }
    }
  }
}
