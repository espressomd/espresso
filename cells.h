/*************************************************/
/*******************  CELLS.H  *******************/
/*************************************************/
#ifndef CELLS_H
#define CELLS_H

#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "grid.h"

/*******************  Structures  *******************/

/** link cell structure. */
typedef struct {
  /* A word about the interacting neighbour cells:\\
     In a 3D lattice each cell has 26 neighbours. Since we deal 
     with pair forces, it is sufficient to calculate only half 
     of the interactions (Newtons law: actio = reactio). For each 
     cell 13 neighbours. This has only to be done for the inner 
     cells. Caution: This implementation needs double sided ghost 
     communication! For single sided ghost communication one 
     would need some ghost-ghost cell interaction as well, which 
     we do not need!   */
  /** number of interacting neighbour cells . */
  int n_neighbours;
  /** interacting neighbour cell list (linear indices) */
  int *neighbours;

  /** number of particles in the cell. */
  int n_particles;
  /** size of particle index array. */
  int max_particles;
  /** particle index array (local index!). */
  int *particles;
} Cell;

/*******************  Variables  *******************/

/** number of linked cells inside the domain of one node (inner cells). */
int n_inner_cells;
/** inner linked cell grid. */
int cell_grid[3];
/** number of linked cells (inner+ghosts). */
int n_cells;
/** linked cell grid with ghost frame. */
int ghost_cell_grid[3];
/** linked cell list. */
Cell *cells;

/** cell size. */
double cell_size[3];
/** inverse cell size. */
double inv_cell_size[3];

/*******************  Functions  *******************/

/** initialize link cell structures. */
void cells_init();

/** sort all particles into cells. */
void sort_particles_into_cells();

/** exit link cell structures. */
void cells_exit();

#endif
