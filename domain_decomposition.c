#include "domain_decomposition.h"

/** Convenient replace for inner cell check. usage: if(IS_INNER_CELL(m,n,o)) {...} */
#define IS_INNER_CELL(m,n,o) \
  ( m > 0 && m < ghost_cell_grid[0] - 1 && \
    n > 0 && n < ghost_cell_grid[1] - 1 && \
    o > 0 && o < ghost_cell_grid[2] - 1 ) 

/** Convenient replace for ghost cell check. usage: if(IS_GHOST_CELL(m,n,o)) {...} */
#define IS_GHOST_CELL(m,n,o) \
  ( m == 0 || m == ghost_cell_grid[0] - 1 || \
    n == 0 || n == ghost_cell_grid[1] - 1 || \
    o == 0 || o == ghost_cell_grid[2] - 1 ) 


/*************************************************/
int pos_to_cell_grid_ind(double pos[3])
{
  int i,cpos[3];
  
  for(i=0;i<3;i++) {
    cpos[i] = (int)((pos[i]-my_left[i])*inv_cell_size[i])+1;

#ifdef PARTIAL_PERIODIC
    if(periodic[i] == 0) {
      if (cpos[i] < 1)                 cpos[i] = 1;
      else if (cpos[i] > cell_grid[i]) cpos[i] = cell_grid[i];
    }
#endif

#ifdef ADDITIONAL_CHECKS
    if(cpos[i] < 1 || cpos[i] >  cell_grid[i]) {
      fprintf(stderr,"%d: illegal cell position cpos[%d]=%d, ghost_grid[%d]=%d for pos[%d]=%f\n",this_node,i,cpos[i],i,ghost_cell_grid[i],i,pos[i]);
      errexit();
    }
#endif

  }
  return get_linear_index(cpos[0],cpos[1],cpos[2], ghost_cell_grid);  
}

/*************************************************/
int pos_to_capped_cell_grid_ind(double pos[3])
{
  int i,cpos[3];
  
  for(i=0;i<3;i++) {
    cpos[i] = (int)((pos[i]-my_left[i])*inv_cell_size[i])+1;

    if (cpos[i] < 1)
      cpos[i] = 1;
    else if (cpos[i] > cell_grid[i])
      cpos[i] = cell_grid[i];
  }
  return get_linear_index(cpos[0],cpos[1],cpos[2], ghost_cell_grid);  
}


/*************************************************/

int max_num_cells_callback(Tcl_Interp *interp, void *_data)
{
  int data = *(int *)_data;
  if (data < 27) {
    Tcl_AppendResult(interp, "WARNING: max_num_cells has to be at least 27. Set max_num_cells = 27!", (char *) NULL);
    data = 27;
  }
  max_num_cells = data;
  mpi_bcast_parameter(FIELD_MAXNUMCELLS);
  mpi_bcast_event(PARAMETER_CHANGED);
  mpi_bcast_event(TOPOLOGY_CHANGED);
  return (TCL_OK);
}
