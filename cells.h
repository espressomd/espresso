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
#define CELL_STRUCTURE_DD 0
#define CELL_STRUCTURE_N2 1
/*@}*/

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
typedef struct Cell {
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

/// List of cells
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

  /** Called when the current cell structure is invalidated because for example the
      box length has changed. This procedure may NOT destroy the old inner and ghost
      cells, but it should free all other organizational data. Note that parameters
      like the box length or the node_grid may already have changed. Therefore
      organizational data has to be stored independently from variables
      that may be changed from outside. */
  void  (*topology_release)();
  /** Initialize the topology. The argument is list of cells, which particles have to be
      sorted into their cells. The particles might not belong to this node.
      This procedure is used when particle data or cell structure has changed and
      the cell structure has to be reinitialized. This also includes setting up the
      cell_structure array. */
  void  (*topology_init)(CellPList *cplist);
  /** Just resort the particles. Used during integration. The particles are stored in
      the cell structure. Domain decomposition can assume for example that particles
      only have to be sent to neighboring nodes. */
  void  (*exchange_and_sort_particles)();
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

/** initialize with a standard cell structure (domain decomposition) */
void cells_init();

/** create a new initialized cell.*/
Cell *alloc_cell();

/** initialize a list of cell pointers */
MDINLINE void init_cellplist(CellPList *cl) {
  cl->n    = 0;
  cl->max  = 0;
  cl->cell = NULL;
}

/** reallocate a list of cell pointers */
MDINLINE void realloc_cellplist(CellPList *cl, int size)
{
  if(size != cl->max) {
    cl->max = size;
    cl->cell = (Cell **) realloc(cl->cell, sizeof(int)*cl->max);
  }
}

/** reinitialize link cell structures. 
 *
 *  It reallocates the cell structure (\ref #cells) and initializes
 *  the contained cell neighbor structure, verlet lists and particle
 *  lists (see \ref init_cell and \ref init_cell_neighbors).
 *
 *  Then it transfers the particles from the old cell structure to the
 *  new one. 
 */
void cells_re_init();

/** called when the topology has changed, so that the cell system can be reinitialized. */
void cells_changed_topology();

/** Calculate and return the total number of particles on this
    node. */
int cells_get_n_particles();

/** debug function to print particle positions: */
void print_particle_positions();

/** debug function to print ghost positions: */
void print_ghost_positions();

/*@}*/

#endif


