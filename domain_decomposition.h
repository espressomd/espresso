#ifndef DOMAIN_DECOMP_H
#define DOMAIN_DECOMP_H
#include "cells.h"

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
/** size of a cell. */
extern double cell_size[3];
/** maximal skin size. */
extern double max_skin;

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

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize the topology. The argument is list of cells, which particles have to be
    sorted into their cells. The particles might not belong to this node.
    This procedure is used when particle data or cell structure has changed and
    the cell structure has to be reinitialized. This also includes setting up the
    cell_structure array. */
void dd_topology_init(CellPList *cl);

/** Called when the current cell structure is invalidated because for example the
    box length has changed. This procedure may NOT destroy the old inner and ghost
    cells, but it should free all other organizational data. Note that parameters
    like the box length or the node_grid may already have changed. Therefore
    organizational data has to be stored independently from variables
    that may be changed from outside. */
void dd_topology_release();

/** Just resort the particles. Used during integration. The particles are stored in
    the cell structure. Domain decomposition can assume for example that particles
    only have to be sent to neighboring nodes. */
void  dd_exchange_and_sort_particles();

/** implements \ref cell_structure::position_to_cell. */
Cell *dd_position_to_cell(double pos[3]);

/** implements \ref cell_structure::position_to_cell. */
Cell *dd_position_to_node(double pos[3]);

/** Callback for setmd maxnumcells (maxnumcells >= 27). 
    see also \ref max_num_cells */
int max_num_cells_callback(Tcl_Interp *interp, void *_data);


/*@}*/

#endif
