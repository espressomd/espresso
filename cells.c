/*************************************************/
/*******************  CELLS.C  *******************/
/*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cells.h"
#include "debug.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "communication.h"

/** increment size of particle buffer */
#define PART_INCREMENT 5
/** half the number of cell neighbours in 3 Dimensions*/
#define MAX_NEIGHBOURS 13

/*******************  Variables  *******************/

/** inner linked cell grid. */
int cell_grid[3];
/** linked cell grid with ghost frame. */
int ghost_cell_grid[3];
/** number of linked cells inside the domain of one node (inner cells). */
int n_inner_cells;

/** number of linked cells (inner+ghosts). */
int n_cells;
/** linked cell list. */
Cell *cells;

/** cell size. */
double cell_size[3];
/** inverse cell size. */
double inv_cell_size[3];

/*******************  privat functions  *******************/

void init_cell_neighbours();
int  is_inner_cell(int i, int adim, int bdim, int cdim);

void cell_memory_info(int verb);
void print_ifield(int *field, int size);
void print_dfield(double *field, int size);

/*******************  exported functions  *******************/

void cells_init() 
{
  int i;
  int node_pos[3];
  CELL_TRACE(fprintf(stderr,"%d: cells_init \n",this_node));
  /* set up dimensions of the cell grid */

  map_node_array(this_node,&node_pos[0],&node_pos[1],&node_pos[2]);     
  for(i=0;i<3;i++) {
    my_left[i]   = node_pos[i]    *local_box_l[i];
    my_right[i]  = (node_pos[i]+1)*local_box_l[i];    
    cell_grid[i] = (int)(local_box_l[i]/max_range);
    if(cell_grid[i] < 1) 
      fprintf(stderr,"%d: cells_init: Less than one cell in direction %d\n",
	      this_node,i);
    
    ghost_cell_grid[i] = cell_grid[i] + 2;
    cell_size[i]       = local_box_l[i]/(double)cell_grid[i];
    inv_cell_size[i]   = (double)cell_grid[i]/local_box_l[i];
  }
  n_inner_cells = cell_grid[0]*cell_grid[1]*cell_grid[2];
  n_cells = ghost_cell_grid[0]*ghost_cell_grid[1]*ghost_cell_grid[2];

  /* there should be a reasonable number of cells only!
     But we will deal with that later... */
  if(this_node==0) {
    if(n_inner_cells > ((max_seen_particle + 1)/nprocs)+1) 
      fprintf(stderr,"0: cells_init: WARNING: More cells per node %d than particles per node %d\n",n_inner_cells,((max_seen_particle + 1)/nprocs)+1);
  }

  /* allocate space for cell structure */
  cells       = (Cell *)malloc(n_cells*sizeof(Cell));
  for(i=0;i<n_cells;i++) {
    cells[i].n_particles=0;
    cells[i].max_particles = PART_INCREMENT;
    cells[i].particles = malloc(cells[i].max_particles*sizeof(int));
    init_cell_neighbours(i);
  }

  /* debuging 
  fflush(stderr);
  MPI_Barrier(MPI_COMM_WORLD) ; 
  for(i=0;i<nprocs;i++) {
    if(i==this_node) { cell_memory_info(0); }
    MPI_Barrier(MPI_COMM_WORLD) ;
  }
  */
}

/*************************************************/

void cells_exit() 
{
  int i;
  CELL_TRACE(if(this_node<2) fprintf(stderr,"%d: cells_exit:\n",this_node));
  for(i=0;i<n_cells;i++) {
    if(cells[i].n_neighbours>0)  free(cells[i].neighbours);
    if(cells[i].max_particles>0) free(cells[i].particles);
  }
  free(cells);
}

/*************************************************/

void sort_particles_into_cells()
{
  
  int i,n;
  int cpos[3], gcg[3];
  int ind;

  CELL_TRACE(fprintf(stderr,"%d: sort_particles_into_cells:\n",this_node));
 
  for(i=0;i<3;i++) gcg[i] = ghost_cell_grid[i];

  /* remove cell particles */
  for(i=0;i<n_cells;i++) cells[i].n_particles=0;
  
  /* particle loop */
  for(n=0;n<n_particles;n++) {
    /* calculate cell index */
    for(i=0;i<3;i++) 
      cpos[i] = (int)((particles[n].p[i]-my_left[i])*inv_cell_size[i])+1;
    ind = get_linear_index(cpos[0],cpos[1],cpos[2],gcg[0],gcg[1],gcg[2]);
    if(ind<0 || ind > n_cells)
      fprintf(stderr,"%d: illigal cell index !(0<%d<%d)\n" ,this_node,ind, n_cells);
    /* Append particle in particle list of that cell */
    if(cells[ind].n_particles >= cells[ind].max_particles-1) 
    realloc_cell_particles(ind,cells[ind].n_particles+(2*PART_INCREMENT));
    cells[ind].particles[cells[ind].n_particles] = n;
    cells[ind].n_particles++;
  }
 
  /* debuging 
  fflush(stderr);
  MPI_Barrier(MPI_COMM_WORLD) ; 
  for(i=0;i<nprocs;i++) {
    if(i==this_node) { 
      fprintf(stderr,"%d: Cell Memory Information:\n",this_node);
      for(i=0;i<n_cells;i++) { 
	if(cells[i].n_particles > 0) 
	  fprintf(stderr,"C%d(np=%d - %d) ",i,cells[i].n_particles,cells[i].max_particles);
      }
      fprintf(stderr,"\n");
      fflush(stderr);
    }
    MPI_Barrier(MPI_COMM_WORLD) ;
  }
  */
}



/*******************  privat functions  *******************/

void init_cell_neighbours(int i)
{
  int j,m,n,o,cnt=0;
  int cg[3],p1[3],p2[3];

  memcpy(cg,ghost_cell_grid,3*sizeof(int));
  if(is_inner_cell(i,cg[0],cg[1],cg[2])) { 
    cells[i].neighbours = malloc(MAX_NEIGHBOURS*sizeof(int));    
    get_grid_pos(i,&p1[0],&p1[1],&p1[2],cg[0],cg[1],cg[2]);
    /* loop through all neighbours */
    for(m=-1;m<2;m++) 
      for(n=-1;n<2;n++) 
	for(o=-1;o<2;o++) {
	  p2[0] = p1[0]+m;   p2[1] = p1[1]+n;   p2[2] = p1[2]+o;
	  j = get_linear_index(p2[0],p2[1],p2[2],cg[0],cg[1],cg[2]);
	  /* take the upper half of all neighbours 
	     and add them to the neighbour list */
	  if(j > i) {
	    cells[i].neighbours[cnt] = j;
	    cnt++;
	  }
	}
    cells[i].n_neighbours = cnt;
  }
  else {
    cells[i].n_neighbours = 0;
    cells[i].neighbours = NULL;
  }   
}

/*************************************************/

int get_linear_index(int a, int b, int c, int adim, int bdim, int cdim)
{
  return (a + adim*(b + bdim*c));
}

/*************************************************/

void get_grid_pos(int i, int *a, int *b, int *c, int adim, int bdim, int cdim)
{
  *a = i % adim;
  i /= adim;
  *b = i % bdim;
  i /= bdim;
  *c = i;
}

/*************************************************/

/** check wether the cell with linear index i is an inner cell or not. */
int  is_inner_cell(int i, int adim, int bdim, int cdim)
{
  int pos[3];
  get_grid_pos(i,&pos[0],&pos[1],&pos[2],adim,bdim,cdim);
  return (pos[0]>0 && pos[0] < adim-1 &&
	  pos[1]>0 && pos[1] < bdim-1 &&
	  pos[2]>0 && pos[2] < cdim-1);
}

/*************************************************/

/** realloc cell particle array.
    Step size for increase and decrease is PART_INCREMENT.
    @param index    linear cell index.
    @param size     desired size of the  particle array.
 */
void realloc_cell_particles(int index, int size)
{
  int incr;

  if( size > cells[index].max_particles) {
    incr = (size-cells[index].max_particles)/PART_INCREMENT;
    cells[index].max_particles += incr*PART_INCREMENT;
  }
  else if( size < (cells[index].max_particles-PART_INCREMENT) ) {
    incr =(size-cells[index].max_particles)/PART_INCREMENT;
    cells[index].max_particles += incr*PART_INCREMENT;
  }
  else incr=0;

  if(incr != 0) {
    cells[index].particles = (int *)
      realloc(cells[index].particles,sizeof(int)*cells[index].max_particles);
  }
}

void cell_memory_info(int verb)
{
  int i;
  
  fprintf(stderr,"%d: Cell Memory Information:\n",this_node);
  fprintf(stderr,"    Inner Cell Grid: "); print_ifield(cell_grid,3);
  fprintf(stderr,"    Ghost Cell Grid: "); print_ifield(ghost_cell_grid,3);
  fprintf(stderr,"    cells: n_cells = %d, n_inner_cells %d\n",
	  n_cells,n_inner_cells); 
  if(verb>0) {
    for(i=0;i<n_cells;i++) {
      fprintf(stderr,"    Cell %d: neighbours: ",i); 
      print_ifield(cells[i].neighbours,cells[i].n_neighbours);
      fprintf(stderr,"    Cell %d: particles: max_particles %d, ",
	      i,cells[i].max_particles);
      print_ifield(cells[i].particles,cells[i].n_particles);
    }
  }
  fflush(stderr);
}

void print_ifield(int *field, int size)
{
  int i;
  fprintf(stderr,"size = %d ",size);
  if(size>0) {
    fprintf(stderr,"{");
    for(i=0;i<size-1;i++) fprintf(stderr,"%d ",field[i]);
    fprintf(stderr,"%d}",field[size-1]);
  }
  fprintf(stderr,"\n");
}

void print_dfield(double *field, int size)
{
  int i;
  fprintf(stderr,"size = %d {",size);
  for(i=0;i<size-1;i++) fprintf(stderr,"%.2f ",field[i]);
  fprintf(stderr,"%.2f}\n",field[size-1]);
}
