/** \file cells.c
 *
 *  This file contains everything related to the link cell
 *  algorithm. This modul strongly interacts with the ghost particle
 *  structure (ghosts.c) and the verlet list algorithm (verlet.c). The
 *  initialization (cells_init()) and cleaning up (cells_exit()) is
 *  done from the integrator (integrate.c).
 *
 *  The domain of a node is split into a 3D cell grid with dimension
 *  cell_grid. Together with one ghost cell layer on each side the
 *  overall dimension of the ghost cell grid is ghost_cell_grid.
 *  You can see a 2D graphical representation of the linked cell grid below. 
 *
 *  \image html linked_cells.gif "Linked cells structure"
 *  \image latex linked_cells.eps "Linked cells structure" \width=6cm
 *
 *  2D representation of a linked cell grid: n_cells = 64, cell_grid =
 *  {4,4}, ghost_cell_grid = {6,6}
 *
 * Each cell has 3^D-1 neighbour cells (For cell 14 they are
 * marked). Since we deal with pair forces, it is sufficient to
 * calculate only half of the interactions (Newtons law: actio =
 * reactio). I have chosen the upper half e.g. all neighbour cells with
 * a higher linear index (For cell 14 they are marked in light
 * blue). Caution: This implementation needs double sided ghost
 * communication! For single sided ghost communication one would need
 * some ghost-ghost cell interaction as well, which we do not need!
 *   */
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
#include "utils.h"

/************************************************
 * defines
 ************************************************/

/** increment size of particle buffer */
#define PART_INCREMENT 5
/** half the number of cell neighbours in 3 Dimensions*/
#define MAX_NEIGHBOURS 13

/************************************************
 * variables
 ************************************************/

int cell_grid[3];
int ghost_cell_grid[3];

/** number of linked cells in nodes spatial domain. */
int n_inner_cells;
int n_cells;
Cell *cells;

/** cell size. 
    Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
double cell_size[3];
/** inverse cell size. */
double inv_cell_size[3];

/************************************************
 * privat functions (headers)
 ************************************************/

/** initializes the interacting neighbour cell list
 *  (cells.neighbours).  the created list of interacting neighbour
 *  cells is used by the verlet algorithm (see verlet.c) to build the
 *  verlet list.
 *
 *  @param i linear index of a cell.  */
void init_cell_neighbours(int i);
 
/** check wether the cell with linear index i is an inner cell or not.
 *  returns 1 if the cell is an inner cell and 0 if it is not.
 *  @param i        linear index of the cell to test.
 *  @param gcg[3]   ghost_cell_grid[3].
 */
int  is_inner_cell(int i, int gcg[3]);

/************************************************
 * functions
 ************************************************/

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
    inv_cell_size[i]   = 1.0 / cell_size[i];
  }
  n_inner_cells = cell_grid[0]*cell_grid[1]*cell_grid[2];
  n_cells = ghost_cell_grid[0]*ghost_cell_grid[1]*ghost_cell_grid[2];

  /* there should be a reasonable number of cells only!
     But we will deal with that later... */
  if(this_node==0) {
    if(n_inner_cells > ((max_seen_particle + 1)/n_nodes)+1) 
      fprintf(stderr,"0: cells_init: WARNING: More cells per node %d than particles per node %d\n",n_inner_cells,((max_seen_particle + 1)/n_nodes)+1);
  }

  /* allocate space for cell structure */
  cells       = (Cell *)malloc(n_cells*sizeof(Cell));
  for(i=0;i<n_cells;i++) {
    cells[i].n_particles=0;
    cells[i].max_particles = PART_INCREMENT;
    cells[i].particles = malloc(cells[i].max_particles*sizeof(int));
    init_cell_neighbours(i);
  }

  CELL_TRACE(fprintf(stderr,"%d: my box: (%.2f, %.2f, %.2f) (%.2f, %.2f, %.2f)\n",
		     this_node,my_left[0],my_left[1],my_left[2],
		     my_right[0],my_right[1],my_right[2]));
  CELL_TRACE(fprintf(stderr,"%d: cell grid (%d+2, %d+2, %d+2) \n",this_node,
		     cell_grid[0],cell_grid[1],cell_grid[2]));
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
  CELL_TRACE(fprintf(stderr,"%d: sort %d paricles \n",this_node,n_particles));
  /* particle loop */
  for(n=0;n<n_particles;n++) {
    /* calculate cell index */
    for(i=0;i<3;i++) 
      cpos[i] = (int)((particles[n].p[i]-my_left[i])*inv_cell_size[i])+1;
    ind = get_linear_index(cpos[0],cpos[1],cpos[2], gcg[0],gcg[1],gcg[2]);
    if(ind<0 || ind > n_cells)
      fprintf(stderr,"%d: illigal cell index !(0<%d<%d)\n" ,this_node,ind, n_cells);
    /* Append particle in particle list of that cell */
    if(cells[ind].n_particles >= cells[ind].max_particles-1) {
      realloc_cell_particles(ind,cells[ind].n_particles + 1);
    }
    cells[ind].particles[cells[ind].n_particles] = n;
    CELL_TRACE(fprintf(stderr,"%d: Part %d at (%.2f, %.2f, %.2f) sorted in cell %d as No. %d\n",this_node,n,particles[n].p[0],particles[n].p[1],particles[n].p[2],ind,cells[ind].n_particles));
    cells[ind].n_particles++;
  }
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


/*******************  privat functions  *******************/

void init_cell_neighbours(int i)
{
  int j,m,n,o,cnt=0;
  int cg[3],p1[3],p2[3];

  memcpy(cg,ghost_cell_grid,3*sizeof(int));
  if(is_inner_cell(i,cg)) { 
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

int  is_inner_cell(int i, int gcg[3])
{
  int pos[3];
  get_grid_pos(i,&pos[0],&pos[1],&pos[2],gcg[0],gcg[1],gcg[2]);
  return (pos[0]>0 && pos[0] < gcg[0]-1 &&
	  pos[1]>0 && pos[1] < gcg[1]-1 &&
	  pos[2]>0 && pos[2] < gcg[2]-1);
}

/*************************************************/

void realloc_cell_particles(int index, int size)
{
  int incr;

  if( size > cells[index].max_particles) {
    incr = (size-cells[index].max_particles)/PART_INCREMENT + 1;
    cells[index].max_particles += incr*PART_INCREMENT;
  }
  else if( size < (cells[index].max_particles-PART_INCREMENT) ) {
    incr =(size-cells[index].max_particles)/PART_INCREMENT;
    cells[index].max_particles += incr*PART_INCREMENT;
  }
  else incr=0;

  if(incr != 0) {
    cells[index].particles = (int *)
      realloc(cells[index].particles, sizeof(int)*cells[index].max_particles);
  }
}
