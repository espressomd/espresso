// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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
#include "ghosts.h"
#include "verlet.h"

/** which cell structure is used */
/*@{*/
#define CELL_STRUCTURE_CURRENT 0
#define CELL_STRUCTURE_DOMDEC  1
#define CELL_STRUCTURE_NSQUARE 2
/*@}*/

/************************************************/
/** \name Data Types */
/************************************************/
/*@{*/

/** A cell is a \ref particleList representing a particle group with
    respect to the integration algorithm.
*/
typedef ParticleList Cell;

/** List of cell pointers. */
typedef struct {
  Cell **cell;
  int n;
  int max;
} CellPList;

/** Describes a cell structure */
typedef struct {
  /** type descriptor */
  int type;

  /** Communicator to exchange ghost cell information. */
  GhostCommunicator ghost_cells_comm;
  /** Communicator to exchange ghost particles. */
  GhostCommunicator exchange_ghosts_comm;
  /** Communicator to update ghost positions. */
  GhostCommunicator update_ghost_pos_comm;
  /** Communicator to collect ghost forces. */
  GhostCommunicator collect_ghost_force_comm;

  ///
  int   (*position_to_node)(double pos[3]);
  ///
  Cell *(*position_to_cell)(double pos[3]);
} CellStructure;

/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** list of all cells. */
extern Cell *cells;
/** size of \ref cells */
extern int n_cells;
/** list of all cells containing particles physically on the local node */
extern CellPList local_cells;
/** list of all cells containing ghosts */
extern CellPList ghost_cells;

/** Type of cell structure in use */
extern CellStructure cell_structure;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize the cell structure on program start with the default cell structure of type domain decomposition. */
void cells_pre_init();

/** Reallocate the list of all cells (\ref cells). */
void realloc_cells(int size);

/** initialize a list of cell pointers */
MDINLINE void init_cellplist(CellPList *cpl) {
  cpl->n    = 0;
  cpl->max  = 0;
  cpl->cell = NULL;
}

/** reallocate a list of cell pointers */
MDINLINE void realloc_cellplist(CellPList *cpl, int size)
{
  if(size != cpl->max) {
    cpl->max = size;
    cpl->cell = (Cell **) realloc(cpl->cell, sizeof(int)*cpl->max);
  }
}

/** reinitialize the cell structures.
    @param new_cs gives the new topology to use afterwards. May be set to
    \ref CELL_STRUCTURE_CURRENT for not changing it.
*/
void cells_re_init(int new_cs);

/** Calculate and return the total number of particles on this
    node. */
int cells_get_n_particles();

/** debug function to print particle positions: */
void print_local_particle_positions();

/** debug function to print ghost positions: */
void print_ghost_positions();

/** implementation of the Tcl command \ref tcl_cellsystem */
int cellsystem(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv);
/*@}*/

#endif


