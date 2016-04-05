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

#include "lattice_inline.hpp"

int lattice_switch = LATTICE_OFF ;

//int Lattice::init(double *agrid, double* offset, int halo_size, size_t dim) {
int Lattice::init(double *agrid, double* offset, int halo_size, size_t dim) {
    this->dim=dim;

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

    this->_data = Utils::malloc(this->element_size*this->halo_grid_volume);
    memset(this->_data, (unsigned int)(-1), this->element_size*this->halo_grid_volume);

}

void Lattice::interpolate(double* pos, double* value) {
    if (this->interpolation_type == INTERPOLATION_LINEAR) {
        interpolate_linear(pos, value);
    } else {
        runtimeErrorMsg() <<"Unknown interpolation type";
    }
}

void Lattice::interpolate_linear(double* pos, double* value) {
    int left_halo_index[3];
    double d[3];
    if (this->halo_size <= 0) {
        runtimeErrorMsg() <<"Error in interpolate_linear: halo size is 0";
        return;
    }
    for (int dim = 0; dim<3; dim++) {
        left_halo_index[dim]=(int) floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]) + this->halo_size;
        d[dim]=((pos[dim]-this->local_offset[dim])/this->agrid[dim] - floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]));
        if (left_halo_index[dim] < 0 || left_halo_index[dim] >= this->halo_grid[dim]) {
            runtimeErrorMsg() <<"Error in interpolate_linear: Particle out of range";
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
        runtimeErrorMsg() <<"Unknown interpolation type";
    }
}

void Lattice::interpolate_linear_gradient(double* pos, double* value) {
    int left_halo_index[3];
    double d[3];
    if (this->halo_size <= 0) {
        runtimeErrorMsg() << "Error in interpolate_linear: halo size is 0";
        return;
    }
    for (int dim = 0; dim<3; dim++) {
        left_halo_index[dim]=(int) floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]) + this->halo_size;
        d[dim]=((pos[dim]-this->local_offset[dim])/this->agrid[dim] - floor((pos[dim]-this->local_offset[dim])/this->agrid[dim]));
        if (left_halo_index[dim] < 0 || left_halo_index[dim] >= this->halo_grid[dim]) {
            runtimeErrorMsg() <<"Error in interpolate_linear: Particle out of range";
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
        global_index[i] = (int)dround((pos[i]-this->offset[i])/this->agrid[i]);
    }

    for (int i=-this->halo_size; i<=this->halo_size; i++) {
        for (int j=-this->halo_size; j<=this->halo_size; j++) {
            for (int k=-this->halo_size; k<=this->halo_size; k++) {
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

int global_pos_in_local_box(double pos[3]) {
    if (!(pos[0]>my_left[0]  &&  pos[0]<my_right[0] ))
        return 0;
    if (!(pos[1]>my_left[1]  &&  pos[1]<my_right[1] ))
        return 0;
    if (!(pos[2]>my_left[2]  &&  pos[2]<my_right[2] ))
        return 0;
    return 1;
}

int Lattice::global_pos_in_local_halobox(double pos[3]) {
    for (unsigned int i=0; i<this->dim; i++) {
        if (!(pos[i]>this->local_offset[i]-this->halo_size*this->agrid[i] &&
              pos[i]<this->local_offset[i]+this->halo_grid[i]*this->agrid[i] ))
            return 0;
    }
    return 1;
}

int Lattice::global_pos_to_lattice_index_checked(double pos[3], int* index) {
    int i;
    for (i=0; i<3; i++)
        if (fabs(fmod(pos[i]-this->offset[i],this->agrid[i])) > ROUND_ERROR_PREC)
            return ES_ERROR;
    int ind[3];
    for (i=0; i<3; i++)
        ind[i] = (int) round((pos[i]-this->offset[i])/this->agrid[i]);
    *index = get_linear_index(this->halo_size + ind[0], this->halo_size + ind[1], this->halo_size + ind[2], this->halo_grid);
    return ES_OK;
}

/* @TODO: Implement! */
int Lattice::map_lattice_to_position(int *ind, int *grid) {
    return 0;
}

void Lattice::map_linear_index_to_global_pos(index_t ind, double* pos) {
    int index_in_halogrid[3];
    get_grid_pos(ind, &index_in_halogrid[0], &index_in_halogrid[1], &index_in_halogrid[2], this->halo_grid);
    pos[0] = this->local_offset[0] + (index_in_halogrid[0] - this->halo_size)*this->agrid[0];
    pos[1] = this->local_offset[1] + (index_in_halogrid[1] - this->halo_size)*this->agrid[1];
    pos[2] = this->local_offset[2] + (index_in_halogrid[2] - this->halo_size)*this->agrid[2];
}

void Lattice::map_local_index_to_pos(index_t* index, double* pos) {
    pos[0] = this->local_offset[0] + (index[0])*this->agrid[0];
    pos[1] = this->local_offset[1] + (index[1])*this->agrid[1];
    pos[2] = this->local_offset[2] + (index[2])*this->agrid[2];
}

int Lattice::map_global_index_to_halo_index(index_t* global_index, index_t* halo_index) {
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

void Lattice::map_halo_index_to_pos(index_t* index_in_halogrid, double* pos) {
    pos[0] = this->local_offset[0] + (index_in_halogrid[0] - this->halo_size)*this->agrid[0];
    pos[1] = this->local_offset[1] + (index_in_halogrid[1] - this->halo_size)*this->agrid[1];
    pos[2] = this->local_offset[2] + (index_in_halogrid[2] - this->halo_size)*this->agrid[2];
}

void Lattice::map_position_to_lattice(const double pos[3], index_t node_index[8], double delta[6]) {

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

void Lattice::get_data_for_halo_index(index_t* ind, void** data) {
    (*data) = ((char*)this->_data) + get_linear_index(ind[0], ind[1], ind[2], this->halo_grid)*this->element_size;
}

void Lattice::get_data_for_linear_index(index_t ind, void** data) {
    (*data) = ((char*)this->_data) + ind*this->element_size;
}

void Lattice::get_data_for_local_index(index_t* ind, void** data) {
    index_t index_in_halogrid[3];
    index_in_halogrid[0] = ind[0]+this->halo_size;
    index_in_halogrid[1] = ind[1]+this->halo_size;
    index_in_halogrid[2] = ind[2]+this->halo_size;
    (*data) = ((char*)this->_data) + get_linear_index(index_in_halogrid[0], index_in_halogrid[1], index_in_halogrid[2], this->halo_grid)*this->element_size;
}

void Lattice::set_data_for_local_halo_grid_index(index_t* ind, void* data) {
    memmove(((char*)this->_data) + get_linear_index(ind[0], ind[1], ind[2],  this->halo_grid)*this->element_size, data, this->element_size);

}

void Lattice::set_data_for_local_grid_index(index_t* ind, void* data) {
    memmove(((char*)this->_data) + get_linear_index(ind[0]+this->halo_size, ind[1]+this->halo_size, ind[2]+this->halo_size,  this->halo_grid)*this->element_size, data, this->element_size);
}

int Lattice::global_pos_to_lattice_halo_index(double* pos, index_t*  ind) {
    for (int i = 0; i<3; i++) {
        ind[i] = (int)dround((pos[i]-this->local_offset[i])/this->agrid[i])+this->halo_size;
        if (ind[i] < 0 || ind[i] >= this->halo_grid[i])
            return 0;
    }
    return 1;
}

/********************** static Functions **********************/

void Lattice::map_position_to_lattice_global (double pos[3], int ind[3], double delta[6], double tmp_agrid) {
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
