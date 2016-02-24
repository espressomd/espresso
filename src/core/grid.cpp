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
/** \file grid.cpp   Domain decomposition for parallel computing.
 *
 *  For more information on the domain decomposition, 
 *  see \ref grid.hpp "grid.h". 
*/
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "utils.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "verlet.hpp"
#include "cells.hpp"
#include "interaction_data.hpp"

/************************************************
 * defines
 ************************************************/    

#define MAX_INTERACTION_RANGE 1e100

/**********************************************
 * variables
 **********************************************/

int node_grid[3] = { 0, 0, 0};
int node_pos[3] = {-1,-1,-1};
#ifdef LEES_EDWARDS
int *node_neighbors;
int *node_neighbor_wrap;
int *node_neighbor_lr;
int  my_neighbor_count;
#else
int node_neighbors[6] = {0, 0, 0, 0, 0, 0};
#endif
int boundary[6]       = {0, 0, 0, 0, 0, 0};
int periodic  = 7;

double box_l[3]       = {1, 1, 1};
double box_l_i[3]     = {1, 1, 1};
double min_box_l;
double local_box_l[3] = {1, 1, 1};
double min_local_box_l;
double my_left[3]     = {0, 0, 0};
double my_right[3]    = {1, 1, 1};

/************************************************************/

void init_node_grid()
{
  grid_changed_n_nodes();
  cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
}

int node_grid_is_set()
{
  return (node_grid[0] > 0);
}

int map_position_node_array(double pos[3])
{
  int i, im[3]={0,0,0};
  double f_pos[3];

  for (i = 0; i < 3; i++)
    f_pos[i] = pos[i];

  fold_position(f_pos, im);

  for (i = 0; i < 3; i++) {
    im[i] = (int)floor(node_grid[i]*f_pos[i]*box_l_i[i]);
    if (im[i] < 0)
      im[i] = 0;
    else if (im[i] >= node_grid[i])
      im[i] = node_grid[i] - 1;
  }

  return map_array_node(im);
}

int calc_node_neighbors(int node)
{

  int  dir, neighbor_count; 
#ifdef LEES_EDWARDS
  node_neighbors     = (int *)Utils::realloc(node_neighbors, 6 * sizeof(int));
  node_neighbor_lr   = (int *)Utils::realloc(node_neighbor_lr, 6 * sizeof(int));
  node_neighbor_wrap = (int *)Utils::realloc(node_neighbor_wrap, 6 * sizeof(int));
#endif  
 
  map_node_array(node,node_pos);
  for(dir=0;dir<3;dir++) {
    int buf;

#ifndef LEES_EDWARDS

    MPI_Cart_shift(comm_cart, dir, -1, &buf, &(node_neighbors[2*dir]));
    MPI_Cart_shift(comm_cart, dir, 1, &buf,  &(node_neighbors[2*dir + 1]));

#else

    /* Writes to node_neighbors[] the integer rank of that neighbor
     * ... the 'buf' stores own rank, which is discarded. */
    if( node_pos[dir] % 2 == 0 ){
        MPI_Cart_shift(comm_cart, dir, -1, &buf, &(node_neighbors[2*dir]));
        MPI_Cart_shift(comm_cart, dir, 1, &buf,  &(node_neighbors[2*dir + 1]));

        node_neighbor_lr[2*dir]   = 0;
        node_neighbor_lr[2*dir+1] = 1;

    }else{
        MPI_Cart_shift(comm_cart, dir, 1, &buf, &(node_neighbors[2*dir]));
        MPI_Cart_shift(comm_cart, dir, -1, &buf,  &(node_neighbors[2*dir + 1]));

        node_neighbor_lr[2*dir]   = 1;
        node_neighbor_lr[2*dir+1] = 0;
    }

#endif

    /* left boundary ? */
    if (node_pos[dir] == 0) {
      boundary[2*dir] = 1;
    }
    else {
      boundary[2*dir] = 0;
    }
    /* right boundary ? */
    if (node_pos[dir] == node_grid[dir]-1) {
      boundary[2*dir+1] = -1;
    }
    else {
      boundary[2*dir+1] = 0;
    }
  }
  
  
  neighbor_count = 6;
#ifdef LEES_EDWARDS
   if( boundary[2] == 1 || boundary[3] == -1 ){

      int x_index, abs_coords[3];

      abs_coords[2] = node_pos[2]; /* z constant */

      if( ( boundary[2] == 1 && boundary[3] != -1) || ( boundary[3] == -1 && boundary[2] != 1) ){
        node_neighbors      = (int *)Utils::realloc(node_neighbors, (neighbor_count + node_grid[0] - 1)*sizeof(int));
        node_neighbor_lr    = (int *)Utils::realloc(node_neighbor_lr, (neighbor_count + node_grid[0] - 1)*sizeof(int));
        node_neighbor_wrap  = (int *)Utils::realloc(node_neighbor_wrap, (neighbor_count + node_grid[0] - 1)*sizeof(int)); 
      }else if( boundary[3] == -1 && boundary[2] == 1){
        node_neighbors      = (int *)Utils::realloc(node_neighbors, (neighbor_count + 2 * node_grid[0] - 2)*sizeof(int));
        node_neighbor_lr    = (int *)Utils::realloc(node_neighbor_lr, (neighbor_count + 2 * node_grid[0] - 2)*sizeof(int));
        node_neighbor_wrap  = (int *)Utils::realloc(node_neighbor_wrap, (neighbor_count + 2 * node_grid[0] - 2)*sizeof(int)); 
      }


      for( x_index = 0; x_index < node_grid[0]; x_index++){
        if( x_index != node_pos[0] ){
            abs_coords[0] = x_index;
            if( x_index > node_pos[0] ){
                if( boundary[2] == 1 ){
                    abs_coords[1] = node_grid[1] - 1; /* wraps to upper y */
                    node_neighbors[neighbor_count] = map_array_node(abs_coords);
                    node_neighbor_lr[neighbor_count]   = 0;
                    neighbor_count++;
                } 
                if( boundary[3] == -1 ){
                    abs_coords[1] = 0;           /* wraps to lower y */
                    node_neighbors[neighbor_count] = map_array_node(abs_coords);
                    node_neighbor_lr[neighbor_count]   =  1;
                    neighbor_count++;
                }
            }else{
                if( boundary[3] == -1 ){
                    abs_coords[1] = 0;           /* wraps to lower y */
                    node_neighbors[neighbor_count] = map_array_node(abs_coords);
                    node_neighbor_lr[neighbor_count]   =  1;
                    neighbor_count++;
                }
                if( boundary[2] == 1 ){
                    abs_coords[1] = node_grid[1] - 1; /* wraps to upper y */
                    node_neighbors[neighbor_count] = map_array_node(abs_coords);
                    node_neighbor_lr[neighbor_count]   = 0;
                    neighbor_count++;
                } 
            }   
        }
      }
   }    
#else
  GRID_TRACE(printf("%d: node_grid %d %d %d, pos %d %d %d, node_neighbors ", this_node, node_grid[0], node_grid[1], node_grid[2], node_pos[0], node_pos[1], node_pos[2]));
#endif


#ifdef LEES_EDWARDS
  my_neighbor_count = neighbor_count;/* set the global neighbor count */
  for(neighbor_count = 0; neighbor_count < my_neighbor_count; neighbor_count++ ){

        if( neighbor_count < 6 ) dir = neighbor_count / 2;
        else dir = 1;
  
        if( boundary[2*dir] == 1 && node_neighbor_lr[neighbor_count] == 0 ){
            node_neighbor_wrap[neighbor_count] =  1;
        }
        else if(boundary[2*dir + 1] == -1 && node_neighbor_lr[neighbor_count] == 1 ){
            node_neighbor_wrap[neighbor_count] = -1;
        }else{
            node_neighbor_wrap[neighbor_count] = 0;
        }

  }
   GRID_TRACE(
   for( dir=0; dir < my_neighbor_count; dir++){ fprintf(stderr, "%d: neighbour %d -->  %d lr: %i wrap: %i\n", this_node, dir,node_neighbors[dir],node_neighbor_lr[dir],node_neighbor_wrap[dir]);}
   );
#endif


  return( neighbor_count );
  
}

void grid_changed_box_l()
{
  int i;

  GRID_TRACE(fprintf(stderr,"%d: grid_changed_box_l:\n",this_node));
  GRID_TRACE(fprintf(stderr,"%d: node_pos %d %d %d\n", this_node, node_pos[0], node_pos[1], node_pos[2]));
  GRID_TRACE(fprintf(stderr,"%d: node_grid %d %d %d\n", this_node, node_grid[0], node_grid[1], node_grid[2]));
  for(i = 0; i < 3; i++) {
    local_box_l[i] = box_l[i]/(double)node_grid[i]; 
    my_left[i]   = node_pos[i]    *local_box_l[i];
    my_right[i]  = (node_pos[i]+1)*local_box_l[i];    
    box_l_i[i] = 1/box_l[i];
  }

  calc_minimal_box_dimensions();

#ifdef GRID_DEBUG
  fprintf(stderr,"%d: local_box_l = (%.3f, %.3f, %.3f)\n",this_node,
	  local_box_l[0],local_box_l[1],local_box_l[2]);
  fprintf(stderr,"%d: coordinates: x in [%.3f, %.3f], y in [%.3f, %.3f], z in [%.3f, %.3f]\n",this_node,
	  my_left[0],my_right[0],my_left[1],my_right[1],my_left[2],my_right[2]);
#endif
}

void grid_changed_n_nodes()
{
  int per[3] = { 1, 1, 1 };
  GRID_TRACE(fprintf(stderr,"%d: grid_changed_n_nodes:\n",this_node));

  MPI_Comm_free(&comm_cart);
  
  MPI_Cart_create(MPI_COMM_WORLD, 3, node_grid, per, 0, &comm_cart);

  MPI_Comm_rank(comm_cart, &this_node);

  MPI_Cart_coords(comm_cart, this_node, 3, node_pos);

  calc_node_neighbors(this_node);

#ifdef GRID_DEBUG
  fprintf(stderr,"%d: node_pos=(%d,%d,%d)\n",this_node,node_pos[0],node_pos[1],node_pos[2]);
  fprintf(stderr,"%d: node_neighbors=(%d,%d,%d,%d,%d,%d)\n",this_node,
	  node_neighbors[0],node_neighbors[1],node_neighbors[2],
	  node_neighbors[3],node_neighbors[4],node_neighbors[5]);
  fprintf(stderr,"%d: boundary=(%d,%d,%d,%d,%d,%d)\n",this_node,
	  boundary[0],boundary[1],boundary[2],boundary[3],boundary[4],boundary[5]);
#endif

  grid_changed_box_l();
}

void calc_minimal_box_dimensions()
{
  int i;
  min_box_l = 2*MAX_INTERACTION_RANGE;
  min_local_box_l = MAX_INTERACTION_RANGE;
  for(i=0;i<3;i++) {
    min_box_l       = dmin(min_box_l, box_l[i]);
    min_local_box_l = dmin(min_local_box_l, local_box_l[i]);
  }
}

void calc_2d_grid(int n, int grid[3])
{
  int i;
  i = (int)sqrt((double)n);
  while(i>=1) {
    if(n%i==0) { grid[0] = n/i; grid[1] = i; grid[2] = 1; return; }
    i--;
  }
}

int map_3don2d_grid(int g3d[3],int g2d[3], int mult[3])
{
  int i,row_dir=-1;
  /* trivial case */
  if(g3d[2]==1) { 
    for(i=0;i<3;i++) mult[i]=1;
    return 2;
  }
  if(g2d[0]%g3d[0] == 0) {
    if(g2d[1]%g3d[1] == 0) {row_dir=2; }
    else if(g2d[1]%g3d[2] == 0) {row_dir=1; g2d[2]=g2d[1]; g2d[1]=1; }
  }
  else if(g2d[0]%g3d[1] == 0) {
    if(g2d[1]%g3d[0]==0) {row_dir=2; i=g2d[0]; g2d[0]=g2d[1]; g2d[1]=i; }
    else if(g2d[1]%g3d[2]==0) {row_dir=0; g2d[2]=g2d[1]; g2d[1]=g2d[0]; g2d[0]=1; }
  }
  else if(g2d[0]%g3d[2] == 0) {
    if(g2d[1]%g3d[0]==0) {row_dir=1; g2d[2]=g2d[0]; g2d[0]=g2d[1]; g2d[1]=1; }
    else if(g2d[1]%g3d[1]==0) {row_dir=0; g2d[2]=g2d[0]; g2d[0]=1; }
  }
  for(i=0;i<3;i++) mult[i]=g2d[i]/g3d[i];
  return row_dir;
}

void rescale_boxl(int dir, double d_new) {
  double scale = (dir-3) ? d_new/box_l[dir] : d_new/box_l[0];
  if (scale < 1.) {
    mpi_rescale_particles(dir,scale);
    if (dir < 3) 
      box_l[dir] = d_new;
    else
      box_l[0] = box_l[1] = box_l[2] = d_new;
    mpi_bcast_parameter(FIELD_BOXL);
  }
  else if (scale > 1.) {
    if (dir < 3) 
      box_l[dir] = d_new;
    else
      box_l[0] = box_l[1] = box_l[2] = d_new;
    mpi_bcast_parameter(FIELD_BOXL);
    mpi_rescale_particles(dir,scale);
  }
}
