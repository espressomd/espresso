#ifndef CELLS_H
#define CELLS_H
/** \file cells.h
 *
 *  This file contains everything related to the link cell
 *  algorithm. 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  This modul strongly interacts with the ghost particle
 *  structure (ghosts.c) and the verlet list algorithm (verlet.c). The
 *  initialization (cells_init()) and cleaning up (cells_exit()) is
 *  called from the integrator (integrate.c).
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
 * Each cell has 3^D neighbor cells (For cell 14 they are
 * marked). Since we deal with pair forces, it is sufficient to
 * calculate only half of the interactions (Newtons law: actio =
 * reactio). I have chosen the upper half e.g. all neighbor cells with
 * a higher linear index (For cell 14 they are marked in light
 * blue). Caution: This implementation needs double sided ghost
 * communication! For single sided ghost communication one would need
 * some ghost-ghost cell interaction as well, which we do not need!
 *
 *  For more information on cells,
 *  see \ref cells.c "cells.c"
 */

#include <tcl.h>
#include "particle_data.h"
#include "verlet.h"

/************************************************/
/** \name Data Types */
/************************************************/
/*@{*/

/** Structure containing information about non bonded interactions
    with particles in a neighbor cell. */
typedef struct {
  /** Just for transparency the index of the neighbor cell. */
  int cell_ind;
  /** Pointer to particle list of neighbor cell. */
  ParticleList *pList;
  /** Verlet list for non bonded interactions of a cell with a neighbor cell. */
  PairList vList;
} IA_Neighbor;

/** Structure containing information of a cell. Contains: cell
    neighbor information, particles in cell.
*/
typedef struct {
  /** number of interacting neighbor cells . 

      A word about the interacting neighbor cells:

      In a 3D lattice each cell has 27 neighbors (including
      itself!). Since we deal with pair forces, it is sufficient to
      calculate only half of the interactions (Newtons law: actio =
      reactio). For each cell 13+1=14 neighbors. This has only to be
      done for the inner cells. 

      Caution: This implementation needs double sided ghost
      communication! For single sided ghost communication one would
      need some ghost-ghost cell interaction as well, which we do not
      need! 

      It follows: inner cells: n_neighbors = 14
      ghost cells:             n_neighbors = 0
  */
  int n_neighbors;
  /** Interacting neighbor cell list  */
  IA_Neighbor *nList;
  /** particle list for particles in the cell. */
  ParticleList pList;
} Cell;

/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** linked cell grid in nodes spatial domain. */
extern int cell_grid[3];
/** linked cell grid with ghost frame. */
extern int ghost_cell_grid[3];

/** number of linked cells (inner+ghosts). */
extern int n_cells;

/** Maximal number of cells per node. In order to avoid memory
 *  problems due to the cell grid one has to specify the maximal
 *  number of \ref #cells . The corresponding callback function is
 *  \ref max_num_cells_callback. If the number of cells \ref n_cells,
 *  defined by \ref ghost_cell_grid is larger than max_num_cells the
 *  cell grid is reduced. max_num_cells has to be larger than 27, e.g
 *  one inner cell.  max_num_cells is initialized with the default
 *  value specified in \ref config.h: \ref CELLS_MAX_NUM_CELLS.
 */
extern int max_num_cells;

/** linked cell list. */
extern Cell *cells;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Pre initialization of the link cell structure. Function called in
    modul initialize.c initialize().  Initializes one cell on each
    node to be able to store the particle data there. */
void cells_pre_init();

/** re initialize link cell structures. 
 *
 *  cells_re_init calculates cell grid parameters: \ref cell_grid,
 *  \ref ghost_cell_grid \ref cell_size, \ref inv_cell_size, \ref
 *  n_cells (via calling \ref calc_cell_grid).
 *
 *  It reallocates the cell structure (\ref #cells) and initializes
 *  the contained cell neighbor structure, verlet lists and particle
 *  lists (see \ref init_cell and \ref init_cell_neighbors).
 *
 *  Then it transfers the particles from the old cell structure to the
 *  new one. 
 */
void cells_re_init();

/** Notify cell code of topology change. Reinits cell cstructure if
  * necesarry (\ref cells_re_init). */
void cells_changed_topology();

/** Calculate and return the total number of particles on this
    node. */
int cells_get_n_particles();

/** Allocate space for a particle.
    \param  id  the identity of the new particle
    \param  pos its position
    \return the new particle structure */
Particle *cells_alloc_particle(int id, double pos[3]);

/** return cell grid index for a position.
    \param pos Position of e.g. a particle.
    \return linear cell grid index. */
int pos_to_cell_grid_ind(double pos[3]);

/** return cell grid index for a position.
    positions out of bounds are capped to the
    nearest valid cell.
    \param pos Position of e.g. a particle.
    \return linear cell grid index. */
int pos_to_capped_cell_grid_ind(double pos[3]);

/** Callback for setmd maxnumcells (maxnumcells >= 27). 
    see also \ref max_num_cells */
int max_num_cells_callback(Tcl_Interp *interp, void *_data);

/** returns true iff cell i is not a ghost cell.
    @param i the cell to test
    @param gcg always ghost_cell_grid
*/
int  is_inner_cell(int i, int gcg[3]);

/** Convenient replace for loops over all particles. */
#define INNER_CELLS_LOOP(m,n,o) \
  for(m=1; m<cell_grid[0]+1; m++) \
    for(n=1; n<cell_grid[1]+1; n++) \
      for(o=1; o<cell_grid[2]+1; o++)

/** Convenient replace for loops over all particles and ghosts. */
#define CELLS_LOOP(m,n,o) \
  for(m=0; m<ghost_cell_grid[0]; m++) \
    for(n=0; n<ghost_cell_grid[1]; n++) \
      for(o=0; o<ghost_cell_grid[2]; o++)

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

/** get the cell index associated with the cell coordinates */
#define CELL_IND(m,n,o) (get_linear_index(m,n,o,ghost_cell_grid))

/** get a pointer to the cell associated with the cell coordinates */
#define CELL_PTR(m,n,o) (&cells[get_linear_index(m,n,o,ghost_cell_grid)])

/** debug function to print particle positions: */
void print_particle_positions();

/** debug function to print ghost positions: */
void print_ghost_positions();

/*@}*/

#endif


