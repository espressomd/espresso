/*******************  CELLS.C  *******************/
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "cells.h"

#define DEBUG

/* increment size of particle buffer */
#define PART_INCREMENT 100

void cells_init() 
{
  int i;

#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d:cells_init \n",this_node);
  fprintf(stderr,"box_lX = %f \n",box_l[0]);
#endif
  

  //for(i=0;i<3;i++) {
  //  cell_grid[i] = (int)(box_l[i]/max_range);
    
  //  if(cell_grid[i] < 1) {
  //    fprintf(stderr,"%d: cells_init: Less than one cell in direction %d\n",this_node,i);
  //  }
    
  // ghost_cell_grid[i] = cell_grid[i] + 2;
  //  cell_size[i]       = box_l[i]/(double)cell_grid[i];
  //  inv_cell_size[i]   = (double)cell_grid[i]/box_l[i];
  //}

  //n_cells = ghost_cell_grid[0]*ghost_cell_grid[1]*ghost_cell_grid[2];

  // cells = (Cell *)malloc(n_cells*sizeof(Cell));

  // for(i=0;i<n_cells;i++) {
  //   cells[i].max_particles = PART_INCREMENT;
  //   *(cells[i].particles) = (int)malloc(PART_INCREMENT*sizeof(int));
  // }

#ifdef DEBUG
  //if(this_node == 0) {
  //  fprintf(stderr,"cell_grid = (%d, %d, %d)\n",
  //	    cell_grid[0],cell_grid[1],cell_grid[2]);
  //  fprintf(stderr,"ghost_cell_grid = (%d, %d, %d) = total %d cells\n",
  //	    ghost_cell_grid[0],ghost_cell_grid[1],ghost_cell_grid[2],n_cells);
  //  fprintf(stderr,"cell_size = (%2e, %2e, %2e)\n",
  //	    cell_size[0],cell_size[1],cell_size[2]);
  //}
#endif
}
