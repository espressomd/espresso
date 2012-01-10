/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file grid.c   Domain decomposition for parallel computing.
 *
 *  For more information on the domain decomposition, 
 *  see \ref grid.h "grid.h". 
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "grid.h"
#include "tcl_interface/grid_tcl.h"
#include "communication.h"
#include "verlet.h"
#include "cells.h"
#include "interaction_data.h"

/************************************************
 * defines
 ************************************************/    

#define MAX_INTERACTION_RANGE 1e100

/**********************************************
 * variables
 **********************************************/

int node_grid[3] = { -1, -1, -1};
int node_pos[3] = {-1,-1,-1};
int node_neighbors[6] = {0, 0, 0, 0, 0, 0};
int boundary[6]   = {0, 0, 0, 0, 0, 0};
int periodic  = 7;

double box_l[3]       = {1, 1, 1};
double box_l_i[3]     = {1, 1, 1};
double min_box_l;
double local_box_l[3] = {1, 1, 1};
double min_local_box_l;
double my_left[3]     = {0, 0, 0};
double my_right[3]    = {1, 1, 1};

/************************************************************/

void setup_node_grid()
{
  if (node_grid[0] < 0) {
    /* auto setup, grid not set */
    calc_3d_grid(n_nodes,node_grid);

    mpi_bcast_parameter(FIELD_NODEGRID);
  }
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

void map_node_array(int node, int pos[3])
{
  get_grid_pos(node, pos, pos + 1, pos + 2, node_grid);
}

int map_array_node(int pos[3]) {
  return get_linear_index(pos[0], pos[1], pos[2], node_grid);
}

void calc_node_neighbors(int node)
{
  int dir,j;
  int n_pos[3];
  
  map_node_array(node,node_pos);
  for(dir=0;dir<3;dir++) {
    for(j=0;j<3;j++) n_pos[j]=node_pos[j];
    /* left neighbor in direction dir */
    n_pos[dir] = node_pos[dir] - 1;
    if(n_pos[dir]<0) n_pos[dir] += node_grid[dir];
    node_neighbors[2*dir]     = map_array_node(n_pos);
    /* right neighbor in direction dir */
    n_pos[dir] = node_pos[dir] + 1;
    if(n_pos[dir]>=node_grid[dir]) n_pos[dir] -= node_grid[dir];
    node_neighbors[(2*dir)+1] = map_array_node(n_pos);
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
}

void grid_changed_box_l()
{
  int i;

  GRID_TRACE(fprintf(stderr,"%d: grid_changed_box_l:\n",this_node));

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
  GRID_TRACE(fprintf(stderr,"%d: grid_changed_n_nodes:\n",this_node));

  calc_node_neighbors(this_node);

#ifdef GRID_DEBUG
  fprintf(stderr,"%d: node_pos=(%d,%d,%d)\n",this_node,node_pos[0],node_pos[1],node_pos[2]);
  fprintf(stderr,"%d: node_neighbors=(%d,%d,%d,%d,%d,%d)\n",this_node,
	  node_neighbors[0],node_neighbors[1],node_neighbors[2],
	  node_neighbors[3],node_neighbors[4],node_neighbors[5]);
  fprintf(stderr,"%d: boundary=(%d,%d,%d,%d,%d,%d)\n",this_node,
	  boundary[0],boundary[1],boundary[2],boundary[3],boundary[4],boundary[5]);
#endif
}

void calc_minimal_box_dimensions()
{
  int i;
  min_box_l = 2*MAX_INTERACTION_RANGE;
  min_local_box_l = MAX_INTERACTION_RANGE;
  for(i=0;i<3;i++) {
    /* #ifdef PARTIAL_PERIODIC  
       if(periodic[i]) {
       min_box_l       = dmin(min_box_l, box_l[i]);
       min_local_box_l = dmin(min_local_box_l, local_box_l[i]);
       }
       #else
    */
    min_box_l       = dmin(min_box_l, box_l[i]);
    min_local_box_l = dmin(min_local_box_l, local_box_l[i]);
    /* #endif */
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

void calc_3d_grid(int n, int grid[3])
{
  int i,j,k,max;
  max = 3*n*n + 1;
  /* generate grid in ascending order */
  for(i=1;i<=n;i++) {
    for(j=i;j<=n;j++) { 
      for(k=j;k<=n;k++) {
	if(i*j*k == n && ((i*i)+(j*j)+(k*k)) < max) {
	  grid[0] = k; grid[1] = j;grid[2] = i;
	  max =  ((i*i)+(j*j)+(k*k));
	}
      }
    }
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
