/** \file cells.c
 *
 *  This file contains everything related to the link cell
 *  algorithm. 
 *
 *  For more information on cells,
 *  see \ref cells.h "cells.h"
 *   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
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
int n_cells;
Cell *cells;
int max_num_cells = 512;

/** number of linked cells in nodes spatial domain. */
int n_inner_cells;
/** cell size. 
    Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
double cell_size[3];
/** inverse cell size. */
double inv_cell_size[3];

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Calculate cell grid dimensions, cell sizes and number of cells.
 *  Calculates the cell grid, based on \ref local_box_l and \ref
 *  max_range. If the number of cells is larger than \ref
 *  max_num_cells, it increases max_range until the number of cells is
 *  smaller or equal \ref max_num_cells. It sets: \ref cell_grid, \ref
 *  ghost_cell_grid, \ref cell_size, \ref inv_cell_size, \ref
 *  n_inner_cells and \ref n_cells.
 */
void calc_cell_grid();

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

/*@}*/
/************************************************************/

void cells_init() 
{
  int i;
  int node_pos[3];
  CELL_TRACE(fprintf(stderr,"%d: cells_init \n",this_node));

  /* set up dimensions of the cell grid */
  calc_cell_grid();
 
  /* this could also go to grid.c ! */
  map_node_array(this_node,node_pos);    
  for(i=0;i<3;i++) {
    my_left[i]   = node_pos[i]    *local_box_l[i];
    my_right[i]  = (node_pos[i]+1)*local_box_l[i];    
  }

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
  if(this_node==0) {
    fprintf(stderr,"0: cell_grid: (%d, %d, %d)\n",cell_grid[0],cell_grid[1],cell_grid[2]);
  }

}


/*************************************************/

void sort_particles_into_cells()
{
  
  int i,n;
  int cpos[3];
  int ind;

  CELL_TRACE(fprintf(stderr,"%d: sort_particles_into_cells:\n",this_node));
 
  /* remove cell particles */
  for(i=0;i<n_cells;i++) cells[i].n_particles=0;
  CELL_TRACE(fprintf(stderr,"%d: sort %d paricles \n",this_node,n_particles));
  /* particle loop */
  for(n=0;n<n_particles;n++) {
    /* calculate cell index */
    for(i=0;i<3;i++) {
      cpos[i] = (int)((particles[n].p[i]-my_left[i])*inv_cell_size[i])+1;
#ifdef PARTIAL_PERIODIC
      if (cpos[i] < 1)
	cpos[i] = 1;
      else if (cpos[i] > cell_grid[i])
	cpos[i] = cell_grid[i];
#endif
    }
    ind = get_linear_index(cpos[0],cpos[1],cpos[2], ghost_cell_grid);
#ifdef ADDITIONAL_CHECKS
    if(ind<0 || ind > n_cells) {
      fprintf(stderr,"%d: illegal cell index %d %d %d, ghost_grid %d %d %d\n", this_node,
	      cpos[0], cpos[1], cpos[2],
	      ghost_cell_grid[0], ghost_cell_grid[1], ghost_cell_grid[2]);
      errexit();
    }
#endif
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
  return (TCL_OK);
}

/*******************  privat functions  *******************/

void calc_cell_grid()
{
  int i;

  /* normal case */
  n_cells=1;
  for(i=0;i<3;i++) {
    ghost_cell_grid[i] = (int)(local_box_l[i]/max_range) + 2;
    n_cells *= ghost_cell_grid[i];
  }

  /* catch case, n_cells > max_num_cells */
  if(n_cells > max_num_cells) {
    int count;
    double cell_range;
    double step;
    double min_box_l;
    double max_box_l;

    min_box_l = dmin(dmin(local_box_l[0],local_box_l[1]),local_box_l[2]);
    max_box_l = dmax(dmax(local_box_l[0],local_box_l[1]),local_box_l[2]);
    step = ((max_box_l/2.0)-max_range)/100; /* Maximal 100 trials! */
    if(step<0.0) {
      fprintf(stderr,"%d: Error: negative step! Something went wrong in calc_cell_grid(). Ask your local Guru\n",this_node);
      errexit();
    }
    cell_range = max_range;

    count = 100;
    while(n_cells > max_num_cells && --count > 0) {
      cell_range += step;
      n_cells=1;
      for(i=0;i<3;i++) {
	ghost_cell_grid[i] = (int)(local_box_l[i]/cell_range) + 2;
	n_cells *= ghost_cell_grid[i];
      }
    }
    if (count == 0) {
      fprintf(stderr, "%d: Error: no suitable cell grid found (max_num_cells was %d)\n",
	      this_node,max_num_cells);
      errexit();
    }
    /* Give information about possible larger skin. */
    cell_range = dmin(min_box_l,cell_range);
    if( cell_range>max_range && this_node==0 ) {
      fprintf(stderr,"Remark: Your parameters would allow a skin of %f instead of your setting %f\n",cell_range-max_cut,skin);
    }
  }

  /* now set all dependent variables */
  n_inner_cells=1;
  for(i=0;i<3;i++) {
    cell_grid[i] = ghost_cell_grid[i]-2;

    /* catch no inner cells case (even though this should never happen!!!) */
    if(cell_grid[i] < 1) {
#ifdef PARTIAL_PERIODIC
      if (!periodic[i])
	cell_grid[i] = 1;
      else
#endif
	{
	  fprintf(stderr,"%d: cells_init: Less than one cell in direction %d\n",
		  this_node,i);
	  errexit();
	}
    }

    n_inner_cells *= cell_grid[i];
    cell_size[i]     = local_box_l[i]/(double)cell_grid[i];
    inv_cell_size[i] = 1.0 / cell_size[i];
  }
}

void init_cell_neighbours(int i)
{
  int j,m,n,o,cnt=0;
  int p1[3],p2[3];

  if(is_inner_cell(i,ghost_cell_grid)) { 
    cells[i].neighbours = malloc(MAX_NEIGHBOURS*sizeof(int));    
    get_grid_pos(i,&p1[0],&p1[1],&p1[2], ghost_cell_grid);
    /* loop through all neighbours */
    for(m = extended[0] ? 0 : -1;
	m < (extended[1] ? 1 :  2); m++) 
      for(n = extended[2] ? 0 : -1;
	  n < (extended[3] ? 1 :  2); n++)
	for(o = extended[4] ? 0 : -1;
	    o < (extended[5] ? 1 :  2); o++) {
	  p2[0] = p1[0]+m;   p2[1] = p1[1]+n;   p2[2] = p1[2]+o;
	  j = get_linear_index(p2[0],p2[1],p2[2], ghost_cell_grid);
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
  get_grid_pos(i,&pos[0],&pos[1],&pos[2],gcg);
  return (pos[0]>0 && pos[0] < gcg[0]-1 &&
	  pos[1]>0 && pos[1] < gcg[1]-1 &&
	  pos[2]>0 && pos[2] < gcg[2]-1);
}

/*************************************************/

