/*************************************************/
/*******************  CELLS.H  *******************/
/*************************************************/
#ifndef CELLS_H
#define CELLS_H

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


/** linked cell grid with ghost frame. */
extern int ghost_cell_grid[3];
/** inner linked cell grid. */
extern int cell_grid[3];

/** number of linked cells (inner+ghosts). */
extern int n_cells;
/** linked cell list. */
extern Cell *cells;

/** cell size. */
extern double cell_size[3];
/** inverse cell size. */
extern double inv_cell_size[3];


/*******************  Functions  *******************/

/** initialize link cell structures. */
void cells_init();

/** sort all particles into inner cells (no ghosts!). */
void sort_particles_into_cells();

/** exit link cell structures. */
void cells_exit();

/** get the linear index from the position (a,b,c) 
    in a 3D grid (adim,bdim,cdim). */
int  get_linear_index(int a, int b, int c, int adim, int bdim, int cdim);

/** get the position (a,b,c) from the linear index
    in a 3D grid (adim,bdim,cdim). */
void get_grid_pos(int i, int *a, int *b, int *c, int adim, int bdim, int cdim);

/** reallocate particle array of cell[index] to size (multiples of PART_INCREMENT). */
void realloc_cell_particles(int index, int size);

#endif
