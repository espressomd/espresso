/*************************************************/
/*******************  CELLS.C  *******************/
/*************************************************/

#include "cells.h"

#define DEBUG

/* increment size of particle buffer */
#define PART_INCREMENT 10

int get_linear_index(int a, int b, int c, int adim, int bdim, int cdim);
void get_grid_pos(int i, int *a, int *b, int *c, int adim, int bdim, int cdim);
void init_cell_neighbours();

void cells_init() 
{
  int i;
  int node_pos[3];
#ifdef DEBUG
  if(this_node < 2) fprintf(stderr,"%d: cells_init \n",this_node);
#endif
  
  map_node_array(this_node,&node_pos[0],&node_pos[1],&node_pos[2]);
  if(this_node < 2) fprintf(stderr,"%d: node_pos (%d,%d,%d)\n",
			    this_node,node_pos[0],node_pos[1],node_pos[2]);
  /* set up dimensions of the cell grid */
  for(i=0;i<3;i++) {
    local_box_l[i] = box_l[i]/((double)processor_grid[i]);
    my_left[i] = node_pos[i]*local_box_l[i];
    my_right[i] = (node_pos[i]+1)*local_box_l[i];    
    cell_grid[i]   = (int)(local_box_l[i]/max_range);
    
    if(cell_grid[i] < 1) {
      fprintf(stderr,"%d: cells_init: Less than one cell in direction %d\n",this_node,i);
    }
        ghost_cell_grid[i] = cell_grid[i] + 2;
    cell_size[i]       = local_box_l[i]/(double)cell_grid[i];
    inv_cell_size[i]   = (double)cell_grid[i]/local_box_l[i];
  }

  n_inner_cells = cell_grid[0]*cell_grid[1]*cell_grid[2];
  n_cells = ghost_cell_grid[0]*ghost_cell_grid[1]*ghost_cell_grid[2];

  if(this_node==0) {
    if(n_inner_cells > (n_total_particles/nprocs)+1) 
      fprintf(stderr,"0: cells_init: WARNING: More cells per node %d than particles per node %d\n",n_inner_cells,(n_total_particles/nprocs)+1);
  }
  cells = (Cell *)malloc(n_cells*sizeof(Cell));
  
  for(i=0;i<n_cells;i++) {
    cells[i].n_particles=0;
    cells[i].max_particles = PART_INCREMENT;
    cells[i].particles = malloc(PART_INCREMENT*sizeof(int));
    cells[i].neighbours = malloc(26*sizeof(int));
  }
  
#ifdef DEBUG
  if(this_node < 2) {
    fprintf(stderr,"cell_grid = (%d, %d, %d)\n",
	    cell_grid[0],cell_grid[1],cell_grid[2]);
    fprintf(stderr,"ghost_cell_grid = (%d, %d, %d) = total %d cells\n",
	    ghost_cell_grid[0],ghost_cell_grid[1],ghost_cell_grid[2],n_cells);
    fprintf(stderr,"cell_size = (%e, %e, %e)\n",
	    cell_size[0],cell_size[1],cell_size[2]);
    fprintf(stderr,"local_box: size (%.2e, %.2e, %.2e)\n",
	    local_box_l[0],local_box_l[1],local_box_l[2]);
    fprintf(stderr,"           from (%.2e, %.2e, %.2e) to (%.2e, %.2e, %.2e)\n",
	    my_left[0],my_left[1],my_left[2],my_right[0],my_right[1],my_right[2]);
  }
#endif
  init_cell_neighbours();

}

void sort_particles_into_cells()
{
  int i,n;
  int cpos[3];
  int ind;
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: sort_particles_into_cells:\n",this_node); 
#endif
  /* particle loop */
  for(n=0;n<n_particles;n++) {
    for(i=0;i<3;i++) 
      cpos[i] = (int)((particles[n].p[i]-my_left[i])*inv_cell_size[i])+1;
    ind = get_linear_index(cpos[0],cpos[1],cpos[2],ghost_cell_grid[0],ghost_cell_grid[1],ghost_cell_grid[2]);
    fprintf(stderr,"%d: Part %d (%.2e, %.2e, %.2e) in cell %d\n",this_node,n,particles[n].p[0],particles[n].p[1],particles[n].p[2],ind);
  }
}

void cells_exit() 
{
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: cells_exit:\n",this_node); 
#endif
}


/*******************  privat functions  *******************/

void init_cell_neighbours()
{
  int i,j,m,n,o;
  int cg[3];
  int p1[3],p2[3];
  int cnt;

  cg[0] = ghost_cell_grid[0];
  cg[1] = ghost_cell_grid[1];
  cg[2] = ghost_cell_grid[2];

  /* loop through all cells */
  for(i=0;i<n_cells;i++) {
    cnt=0;
    get_grid_pos(i,&p1[0],&p1[1],&p1[2],cg[0],cg[1],cg[2]);
    /* is it an inner cell ? */
    if(p1[0]>0 && p1[0] < cg[0]-1 &&
       p1[1]>0 && p1[1] < cg[1]-1 &&
       p1[2]>0 && p1[2] < cg[2]-1) {
#ifdef DEBUG
      if(this_node==0) 
	fprintf(stderr,"Cell %d pos(%d,%d,%d) with neighbours\n",
		i,p1[0],p1[1],p1[2]);
#endif
      for(m=-1;m<2;m++) 
	for(n=-1;n<2;n++) 
	  for(o=-1;o<2;o++) {
	    p2[0] = p1[0]+m;
	    p2[1] = p1[1]+n;
	    p2[2] = p1[2]+o;
	    j = get_linear_index(p2[0],p2[1],p2[2],cg[0],cg[1],cg[2]);
	    if(j != i) {
	      cells[i].particles[cnt] = j;
#ifdef DEBUG
	      if(this_node==0) fprintf(stderr,"%d ",cells[i].particles[cnt]);
#endif
	      cnt++;
	    }
	  }
#ifdef DEBUG
       if(this_node==0) fprintf(stderr,"\n");     
#endif
    }
    cells[i].n_neighbours = cnt; 
  }
}


int get_linear_index(int a, int b, int c, int adim, int bdim, int cdim)
{
  return (a + adim*(b + bdim*c));
}

void get_grid_pos(int i, int *a, int *b, int *c, int adim, int bdim, int cdim)
{
  *a = i % adim;
  i /= adim;
  *b = i % bdim;
  i /= bdim;
  *c = i;
}



