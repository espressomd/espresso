/************************************************/
/*******************  GRID.C  *******************/
/************************************************/
#include "grid.h"
#include "debug.h"

/**********************************************
 * local settings
 **********************************************/

int processor_grid[3] = { -1, -1, -1};
int pe_pos[3] = {-1,-1,-1};
int neighbors[6] = {0, 0, 0, 0, 0, 0};

/**********************************************
 * procedures
 **********************************************/

/* callback for procgrid */
int pgrid_callback(Tcl_Interp *interp, void *_data)
{
  int *data = (int *)_data;
  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }
  processor_grid[0] = data[0];
  processor_grid[1] = data[1];
  processor_grid[2] = data[2];
  if (!setup_processor_grid()) {
    Tcl_AppendResult(interp, "processor grid does not fit nprocs",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  changed_topology();

  return (TCL_OK);
}

int setup_processor_grid()
{
  if (processor_grid[0] < 0) {
    fprintf(stderr, "not implemented: setup_processor_grid()\n");
    processor_grid[0] = nprocs;
    processor_grid[1] = processor_grid[2] = 1;
    return 1;
  }
  if (processor_grid[0]*processor_grid[1]*processor_grid[2] != nprocs)
    return 0;
  return 1;
}

int processor_grid_is_set()
{
  return (processor_grid[1] > 0);
}

int find_node(double pos[3])
{
  int i;
  double f_pos[3];
  
  for(i=0;i<3;i++) f_pos[i] = pos[i] - floor(pos[i]/box_l[i])*box_l[i];
  GRID_TRACE(fprintf(stderr,"(%e, %e, %e) folded to (%e,%e, %e)\n",
		     pos[0],pos[1],pos[2],f_pos[0],f_pos[1],f_pos[2]));
  return map_array_node((int)floor(processor_grid[0]*f_pos[0]/box_l[0]),
			(int)floor(processor_grid[1]*f_pos[1]/box_l[1]),
			(int)floor(processor_grid[2]*f_pos[2]/box_l[2]));
}

void map_node_array(int node, int *a, int *b, int *c)
{
  *a = node % processor_grid[0];
  node /= processor_grid[0];
  *b = node % processor_grid[1];
  node /= processor_grid[1];
  *c = node;
}

int map_array_node(int a, int b, int c) {
  return (a + processor_grid[0]*(b + processor_grid[1]*c));
}

void calc_neighbors(int node)
{
  int dir,j;
  int n_pos[3];
  
  map_node_array(node,&pe_pos[0],&pe_pos[1],&pe_pos[2]);

  for(dir=0;dir<3;dir++) {
    for(j=0;j<3;j++) n_pos[j]=pe_pos[j];
    /* left neighbor in direction dir */
    n_pos[dir] -= 1;
    if(n_pos[dir]<0) n_pos[dir] += processor_grid[dir];
    neighbors[2*dir]     = map_array_node(n_pos[0],n_pos[1],n_pos[2]);
    /* right neighbor in direction dir */
    n_pos[dir] += 2;
    if(n_pos[dir]>=processor_grid[dir]) n_pos[dir] -= processor_grid[dir];
    neighbors[(2*dir)+1] = map_array_node(n_pos[0],n_pos[1],n_pos[2]);
  }
}
