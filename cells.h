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

/************************************************
 * defines
 ************************************************/    

#define CELLS_FLAG_START   -1
#define CELLS_FLAG_PRE_INIT 0
#define CELLS_FLAG_RE_INIT  1

/************************************************
 * data types
 ************************************************/

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

      It follows: inner cells: n_neighbors=14
      ghost cells: n_neighbors=0
  */
  int n_neighbors;
  /** Interacting neighbor cell list  */
  IA_Neighbor *nList;
  /** particle list for particles in the cell. */
  ParticleList pList;
} Cell;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** linked cell grid in nodes spatial domain. */
extern int cell_grid[3];
/** linked cell grid with ghost frame. */
extern int ghost_cell_grid[3];

/** number of linked cells (inner+ghosts). */
extern int n_cells;
/** linked cell list. */
extern Cell *cells;

/** Maximal number of cells per node. 
 *  In order to avoid memory problems due to the cell grid one has to
 *  specify the maximal number of \ref cells. The corresponding
 *  callback function is \ref max_num_cells_callback. If the number of
 *  cells \ref n_cells, defined by \ref ghost_cell_grid is larger than
 *  max_num_cells the cell grid is reduced. max_num_cells has to be
 *  larger than 27, e.g one inner cell.
 */
extern int max_num_cells;

/** cell initialization status. 
    initialized:      cells_init_flag = CELLS_FLAG_START.
    cells_pre_init(): cells_init_flag = CELLS_FLAG_PRE_INIT
    cells_init():     cells_init_flag = CELLS_FLAG_RE_INIT
    cells_exit():     ruft am ende cells_pre_init() auf.
 */
extern int cells_init_flag;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/


/** Pre initialization of the link cell structure. 

Function called in modul initialize.c initialize().  Initializes one
cell on each node to be able to store the particle data there.
*/
void cells_pre_init();

/** Notify cell code of topology change. 

Recalculates the cell sizes.
*/
void cells_changed_topology();

/** re initialize link cell structures. 
 *
 *  cells_init calculates the cell grid dimension (cell_grid[3],
 *  ghost_cell_grid[3]) and allocates the space for the cell
 *  structure: cells.  Further it calculates the values for:
 *  cell_size[3], inv_cell_size[3], n_cells, n_inner_cells.
 *
 *  At the moment it also calculates the edges of the nodes
 *  spatial domain: my_left[3] and my_right[3].
 *
 *  Then it allocates space for the particle index list of each cell
 *  (cells[i].particles) with size PART_INCREMENT and initializes the
 *  neighbor list for the cells (init_cell_neighbors()).  
 */
void cells_re_init();

/** sort all particles into inner cells (no ghosts!). 
 *
 *  In order to build the verlet list (verlet.c) from the link cell
 *  structure one has to sort the particles into the cells. This is
 *  done after the particle exchange (exchange_part()).  
 *
 *  Sorting: Go through local particle list (0...n_particles) and
 *  store the local index into the particle list of the corresponding
 *  cell. 
 *
 *  cell particle index buffer (cells[i].particles) is reallocated if
 *  necessary. ATTENTION: there is at the moment no routine that
 *  reduces the size of this list if it would be possible. You have to
 *  use cells_exit() for that and than reinitialize the cell structure
 *  with cells_init() again. */
void sort_particles_into_cells();

/** exit link cell structures.  
    free space for linked cell structure.  */
void cells_exit();

/** Callback for setmd maxnumcells (maxnumcells >= 27). 
    see also \ref max_num_cells */
int max_num_cells_callback(Tcl_Interp *interp, void *_data);

#endif

/*@}*/
