/*************************************************/
/*******************  CELLS.H  *******************/
/*************************************************/
#ifndef CELLS_H
#define CELLS_H
#include "global.h"


/*******************  Structures  *******************/

/** link cell structure. */
typedef struct {
  /** Position inside the cell grid. */
  double pos[3];
  
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

/** inner linked cell grid. */
int cell_grid[3];
/** linked cell grid with ghost frame. */
int ghost_cell_grid[3];
/** number of linked cells. */
int n_cells;
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


#endif
