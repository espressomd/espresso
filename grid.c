/** \file grid.c
 *
 *  Domain decomposition for parallel computing.
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different nodes (or processes
 *  or processors or threads if you want). This grid is described in
 *  \ref node_grid. Each node has a number \ref this_node and a
 *  position \ref node_pos in that grid. Each node has also 6 nearest
 *  neighbors \ref node_neighbors which are necessary for the
 *  communication between the nodes (see also \ref ghosts.c and \ref
 *  p3m.c for more details about the communication.
 *
 *  For the 6 directions \anchor directions we have the following convention:
 *
 *  \image html directions.gif "Convention for the order of the directions"
 *  \image latex directions.eps "Convention for the order of the directions" width=6cm
 *
 *  The Figure illustrates the direction convetion used for arrays
 *  with 6 (e.g. \ref node_neighbors, \ref boundary) and 3 entries
 *  (e.g \ref node_grid, \ref box_l , \ref my_left,...).
 *  
 *  Attention: If you change anything of the simulation box dimensions
 *  you have to call \ref changed_topology.
 * */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "grid.h"
#include "utils.h"
#include "debug.h"
#include "communication.h"
#include "verlet.h"

/* Rebuild the verlet list when topology changes. */

/**********************************************
 * variables
 **********************************************/

int node_grid[3] = { -1, -1, -1};
int node_pos[3] = {-1,-1,-1};
int node_neighbors[6] = {0, 0, 0, 0, 0, 0};
int boundary[6]   = {0, 0, 0, 0, 0, 0};
int extended[6] = {0, 0, 0, 0, 0, 0};
int periodic[3]    = {1, 1, 1};

double box_l[3]       = {1, 1, 1};
double local_box_l[3] = {1, 1, 1};
double my_left[3]     = {0, 0, 0};
double my_right[3]    = {1, 1, 1};

/**********************************************
 * procedures
 **********************************************/

/* callback for box_l change */
int boxl_callback(Tcl_Interp *interp, void *_data)
{
  double *data = _data;

  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }

  box_l[0] = data[0];
  box_l[1] = data[1];
  box_l[2] = data[2];

  changed_topology();

  return (TCL_OK);
}

#ifdef PARTIAL_PERIODIC
/* callback for periodic change */
int per_callback(Tcl_Interp *interp, void *_data)
{
  int i;
  int *data = _data;

  for (i = 0; i < 3; i++)
    periodic[i] = (data[i] != 0);

  mpi_bcast_parameter(FIELD_PERIODIC);

  return (TCL_OK);
}
#endif

/* callback for node grid */
int node_grid_callback(Tcl_Interp *interp, void *_data)
{
  int *data = (int *)_data;
  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }

  if (data[0]*data[1]*data[2] != n_nodes) {
    Tcl_AppendResult(interp, "node grid does not fit n_nodes",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  node_grid[0] = data[0];
  node_grid[1] = data[1];
  node_grid[2] = data[2];

  changed_topology();

  return (TCL_OK);
}

void setup_node_grid()
{
  if (node_grid[0] < 0) {
    /* auto setup, grid not set */
    fprintf(stderr, "not implemented: setup_node_grid()\n");
    node_grid[0] = n_nodes;
    node_grid[1] = node_grid[2] = 1;

    changed_topology();
  }
}

int node_grid_is_set()
{
  return (node_grid[1] > 0);
}

void changed_topology()
{
  int i;
  for(i = 0; i < 3; i++) {
    local_box_l[i] = box_l[i]/(double)node_grid[i]; 
  }

  rebuild_verletlist = 1;

  /* a little bit of overkill, but safe */
  mpi_bcast_parameter(FIELD_NGRID);
  mpi_bcast_parameter(FIELD_BOXL);
  mpi_bcast_parameter(FIELD_LBOXL);
  mpi_bcast_parameter(FIELD_VERLET);
}

int find_node(double pos[3])
{
  int i, im[3];
  double f_pos[3];

  for (i = 0; i < 3; i++)
    f_pos[i] = pos[i];
  
  fold_particle(f_pos, im);

  for (i = 0; i < 3; i++) {
    im[i] = (int)floor(node_grid[i]*f_pos[i]/box_l[i]);
#ifdef PARTIAL_PERIODIC
    if (!periodic[i]) {
      if (im[i] < 0)
	im[i] = 0;
      else if (im[i] >= node_grid[i])
	im[i] = node_grid[i] - 1;
    }
#endif
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
      extended[2*dir] = periodic[dir] ? 0 : 1;
    }
    else {
      boundary[2*dir] =
	extended[2*dir] = 0;
    }
    /* right boundary ? */
    if (node_pos[dir] == node_grid[dir]-1) {
      boundary[2*dir+1] = -1;
      extended[2*dir+1] = periodic[dir] ? 0 : 1;
    }
    else {
      boundary[2*dir+1] =
	extended[2*dir+1] = 0;
    }
  }
}
