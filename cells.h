/*************************************************/
/*******************  CELLS.H  *******************/
/*************************************************/
#ifndef CELLS_H
#define CELLS_H

/*******************  Structures  *******************/

/** link cell structure. */
typedef struct {
  /** number of interacting neighbour cells . 
      A word about the interacting neighbour cells:\\
      In a 3D lattice each cell has 26 neighbours. Since we deal 
      with pair forces, it is sufficient to calculate only half 
      of the interactions (Newtons law: actio = reactio). For each 
      cell 13 neighbours. This has only to be done for the inner 
      cells. Caution: This implementation needs double sided ghost 
      communication! For single sided ghost communication one 
      would need some ghost-ghost cell interaction as well, which 
      we do not need! */
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

/** initialize link cell structures. 

    cells_init calculates the cell grid dimension and allocates the
    space for the cell structure: cells. 

    Then it allocates space for the particle index list of each cell
    (cells[i].particles) and initializes the neighbour list for the
    cells (init_cell_neighbours()).  */
void cells_init();

/** sort all particles into inner cells (no ghosts!). 

    In order to build the verlet list (verlet.c) from the link cell
    structure one has to sort the particles into the cells. This is
    done after the particle exchange (exchange_part()).  

    Sorting: Go through local particle list (0...n_particles) and
    store the local index into the particle list of the corresponding
    cell. */
void sort_particles_into_cells();

/** exit link cell structures.  
    free space for linked cell structure.  */
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
