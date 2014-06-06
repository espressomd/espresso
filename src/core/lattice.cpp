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
/** \file lattice.cpp 
 *
 * Lattice data structures
 *
 */

#include "utils.hpp"
#include "grid.hpp"
#include "lattice.hpp"

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
void _lattice_allocate_memory(Lattice *lattice); 
int init_lattice(Lattice *lattice, double *agrid, double* offset, int halo_size, size_t dim) {
  int dir;

  lattice->dim=dim;

 /* determine the number of local lattice nodes */
  for (int d=0; d<3; d++) {
    lattice->agrid[d] = agrid[d];
    lattice->global_grid[d] = (int)floor(box_l[d]/agrid[d]+ROUND_ERROR_PREC);
    lattice->offset[d]=offset[d];
    lattice->local_index_offset[d]=(int) ceil((my_left[d]-lattice->offset[d])/lattice->agrid[d]);
    lattice->local_offset[d] = lattice->offset[d] +
      lattice->local_index_offset[d]*lattice->agrid[d];
    lattice->grid[d] = (int) ceil ( ( my_right[d] - lattice->local_offset[d]-ROUND_ERROR_PREC )
      / lattice->agrid[d]);
  }

  /* sanity checks */
  for (dir=0;dir<3;dir++) {
    /* check if local_box_l is compatible with lattice spacing */
//    if (fabs(local_box_l[dir]-lattice->grid[dir]*agrid[dir]) > ROUND_ERROR_PREC*box_l[dir]) {
//      char *errtxt = runtime_error(128);
//      ERROR_SPRINTF(errtxt, "{097 Lattice spacing agrid[%d]=%f is incompatible with local_box_l[%d]=%f (box_l[%d]=%f node_grid[%d]=%d) %f} ",dir,agrid[dir],dir,local_box_l[dir],dir,box_l[dir],dir,node_grid[dir],local_box_l[dir]-lattice->grid[dir]*agrid[dir]);
//      return ES_ERROR;
//    }
  /* set the lattice spacing */
  }

  lattice->element_size = lattice->dim*sizeof(double);

  LATTICE_TRACE(fprintf(stderr,"%d: box_l (%.3f,%.3f,%.3f) grid (%d,%d,%d) node_neighbors (%d,%d,%d,%d,%d,%d)\n",this_node,local_box_l[0],local_box_l[1],local_box_l[2],lattice->grid[0],lattice->grid[1],lattice->grid[2],node_neighbors[0],node_neighbors[1],node_neighbors[2],node_neighbors[3],node_neighbors[4],node_neighbors[5]));

  lattice->halo_size = halo_size;
  /* determine the number of total nodes including halo */
  lattice->halo_grid[0] = lattice->grid[0] + 2*halo_size ;
  lattice->halo_grid[1] = lattice->grid[1] + 2*halo_size ;
  lattice->halo_grid[2] = lattice->grid[2] + 2*halo_size ;

  lattice->grid_volume = lattice->grid[0]*lattice->grid[1]*lattice->grid[2] ;
  lattice->halo_grid_volume = lattice->halo_grid[0]*lattice->halo_grid[1]*lattice->halo_grid[2] ;
  lattice->halo_grid_surface = lattice->halo_grid_volume - lattice->grid_volume ;
  lattice->halo_offset = get_linear_index(halo_size,halo_size,halo_size,lattice->halo_grid) ;
  
  lattice->interpolation_type = INTERPOLATION_LINEAR;
  
  _lattice_allocate_memory(lattice);
  return ES_OK;

}

void _lattice_allocate_memory(Lattice *lattice) {

  lattice->_data = malloc(lattice->element_size*lattice->halo_grid_volume);
  memset(lattice->_data, (unsigned int)(-1), lattice->element_size*lattice->halo_grid_volume);

}

void lattice_interpolate_linear(Lattice* lattice, double* pos, double* value); 

void lattice_interpolate(Lattice* lattice, double* pos, double* value) {
  if (lattice->interpolation_type == INTERPOLATION_LINEAR) {
    lattice_interpolate_linear(lattice, pos, value);
  } else {
    char* c = runtime_error(128);
    ERROR_SPRINTF(c, "Unknown interpolation type");
  }
}

void lattice_interpolate_linear_gradient(Lattice* lattice, double* pos, double* value);

void lattice_interpolate_gradient(Lattice* lattice, double* pos, double* value) {
  if (lattice->interpolation_type == INTERPOLATION_LINEAR) {
    lattice_interpolate_linear_gradient(lattice, pos, value);
  } else {
    char* c = runtime_error(128);
    ERROR_SPRINTF(c, "Unknown interpolation type");
  }
}


void lattice_interpolate_linear_gradient(Lattice* lattice, double* pos, double* value) {
  int left_halo_index[3];
  double d[3];
   if (lattice->halo_size <= 0) {
     char* c = runtime_error(128);
     ERROR_SPRINTF(c, "Error in lattice_interpolate_linear: halo size is 0");
     return;
   }
   for (int dim = 0; dim<3; dim++) {
     left_halo_index[dim]=(int) floor((pos[dim]-lattice->local_offset[dim])/lattice->agrid[dim]) + lattice->halo_size;
     d[dim]=((pos[dim]-lattice->local_offset[dim])/lattice->agrid[dim] - floor((pos[dim]-lattice->local_offset[dim])/lattice->agrid[dim]));
     if (left_halo_index[dim] < 0 || left_halo_index[dim] >= lattice->halo_grid[dim]) {
       char* c = runtime_error(128);
       ERROR_SPRINTF(c, "Error in lattice_interpolate_linear: Particle out of range");
       return;
     }
   }
   
   index_t index;
   double* local_value;

   for (unsigned int i = 0; i<3*lattice->dim; i++) {
     value[i] = 0;
   }
   
   index=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2], lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*(1-d[1])*(1-d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= (1-d[0])*( -1   )*(1-d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= (1-d[0])*(1-d[1])*(  -1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2], lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*(1-d[1])*(1-d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= ( +d[0])*( -1   )*(1-d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= ( +d[0])*(1-d[1])*(  -1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2], lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*( +d[1])*(1-d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= (1-d[0])*( +1   )*(1-d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= (1-d[0])*( +d[1])*(  -1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2], lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*( +d[1])*(1-d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= ( +d[0])*( +1   )*(1-d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= ( +d[0])*( +d[1])*(  -1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]  , left_halo_index[1]  , left_halo_index[2] + 1, lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*(1-d[1])*( +d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= (1-d[0])*( -1   )*( +d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= (1-d[0])*(1-d[1])*(  +1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2]+1, lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*(1-d[1])*( +d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= ( +d[0])*( -1   )*( +d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= ( +d[0])*(1-d[1])*(  +1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2]+1, lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  -1  )*( +d[1])*( +d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= (1-d[0])*( +1   )*( +d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= (1-d[0])*( +d[1])*(  +1  ) * local_value[i] / lattice->agrid[2];
   }
   index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2]+1, lattice->halo_grid);
   for (unsigned int i = 0; i<lattice->dim; i++) {
     lattice_get_data_for_linear_index(lattice, index, (void**) &local_value);
     value[3*i  ]+= (  +1  )*( +d[1])*( +d[2]) * local_value[i] / lattice->agrid[0];
     value[3*i+1]+= ( +d[0])*( +1   )*( +d[2]) * local_value[i] / lattice->agrid[1];
     value[3*i+2]+= ( +d[0])*( +d[1])*(  +1  ) * local_value[i] / lattice->agrid[2];
   }

}


void lattice_interpolate_linear(Lattice* lattice, double* pos, double* value) {
  int left_halo_index[3];
  double d[3];
   if (lattice->halo_size <= 0) {
     char* c = runtime_error(128);
     ERROR_SPRINTF(c, "Error in lattice_interpolate_linear: halo size is 0");
     return;
   }
   for (int dim = 0; dim<3; dim++) {
     left_halo_index[dim]=(int) floor((pos[dim]-lattice->local_offset[dim])/lattice->agrid[dim]) + lattice->halo_size;
     d[dim]=((pos[dim]-lattice->local_offset[dim])/lattice->agrid[dim] - floor((pos[dim]-lattice->local_offset[dim])/lattice->agrid[dim]));
     if (left_halo_index[dim] < 0 || left_halo_index[dim] >= lattice->halo_grid[dim]) {
       char* c = runtime_error(128);
       ERROR_SPRINTF(c, "Error in lattice_interpolate_linear: Particle out of range");
       return;
     }
   }
   double w[8];
   index_t index[8];
   w[0] = (1-d[0])*(1-d[1])*(1-d[2]);
   index[0]=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2], lattice->halo_grid);
   w[1] = ( +d[0])*(1-d[1])*(1-d[2]);
   index[1]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2], lattice->halo_grid);
   w[2] = (1-d[0])*( +d[1])*(1-d[2]);
   index[2]=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2], lattice->halo_grid);
   w[3] = ( +d[0])*( +d[1])*(1-d[2]);
   index[3]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2], lattice->halo_grid);

   w[4] = (1-d[0])*(1-d[1])*( +d[2]);
   index[4]=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2]+1, lattice->halo_grid);
   w[5] = ( +d[0])*(1-d[1])*( +d[2]);
   index[5]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2]+1, lattice->halo_grid);
   w[6] = (1-d[0])*( +d[1])*( +d[2]);
   index[6]=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2]+1, lattice->halo_grid);
   w[7] = ( +d[0])*( +d[1])*( +d[2]);
   index[7]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2]+1, lattice->halo_grid);
   
   for (unsigned int i = 0; i<lattice->dim; i++) {
     value[i] = 0;
   }

   double* local_value;
   for (unsigned int i=0; i<8; i++) {
     lattice_get_data_for_linear_index(lattice, index[i], (void**) &local_value);
     for (unsigned int j = 0; j<lattice->dim; j++) {
       value[j]+=w[i]*local_value[j];
     }
   }
}

void lattice_set_data_for_global_position_with_periodic_image(Lattice* lattice, double* pos, void* data) {

  index_t replica[3];
  index_t global_index[3];
  index_t halo_index[3];

  
  for (int i = 0; i<3; i++) {
    global_index[i] = (int) floor((pos[i]-lattice->offset[i])/lattice->agrid[i]+ROUND_ERROR_PREC);
  }

  for (int i=-2; i<=2; i++) {
    for (int j=-2; j<=2; j++) {
      for (int k=-2; k<=2; k++) {
        replica[0]=global_index[0]+i*lattice->global_grid[0];
        replica[1]=global_index[1]+j*lattice->global_grid[1];
        replica[2]=global_index[2]+k*lattice->global_grid[2];
        if (map_global_index_to_halo_index( lattice, replica, halo_index)) {
          lattice_set_data_for_local_halo_grid_index(lattice, halo_index, data);
        }

      }
    }
  }
}


