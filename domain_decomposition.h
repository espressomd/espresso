// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef DOMAIN_DECOMP_H
#define DOMAIN_DECOMP_H

/** \file domain_decomposition.h
 *
 *  This file contains everything related to the cell system: domain decomposition.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  The simulation box is split into spatial domains for each node
 *  according to a cartesian node grid (\ref node_grid).
 *
 *  The domain of a node is split into a 3D cell grid with dimension
 *  \ref DomainDecomposition::cell_grid. Together with one ghost cell
 *  layer on each side the overall dimension of the ghost cell grid is
 *  \ref DomainDecomposition::ghost_cell_grid. The domain
 *  decomposition enables one the use of the linked cell algorithm
 *  which is in turn used for setting up the verlet list for the
 *  system. You can see a 2D graphical representation of the linked
 *  cell grid below.
 *
 *  \image html  linked_cells.gif "Linked cells structure"
 *  \image latex linked_cells.eps "Linked cells structure" 
 *
 *  2D representation of a linked cell grid: cell_grid =
 *  {4,4}, ghost_cell_grid = {6,6}
 *
 * Each cell has 3^D neighbor cells (For cell 14 they are
 * marked). Since we deal with pair forces, it is sufficient to
 * calculate only half of the interactions (Newtons law: actio =
 * reactio). We have chosen the upper half e.g. all neighbor cells with
 * a higher linear index (For cell 14 they are marked in light
 * blue). Caution: This implementation needs double sided ghost
 * communication! For single sided ghost communication one would need
 * some ghost-ghost cell interaction as well, which we do not need!
 *
 *  For more information on cells,
 *  see \ref cells.h 
*/

#include "cells.h"
#include "integrate.h"
#include "communication.h"
#include "verlet.h"
#include "thermostat.h"
#include "debug.h"

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
} IA_Neighbor_List;

/** Structure containing the information about the cell grid used for domain decomposition. */
typedef struct {
  /** linked cell grid in nodes spatial domain. */
  int cell_grid[3];
  /** linked cell grid with ghost frame. */
  int ghost_cell_grid[3];
  /** cell size. 
      Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
  double cell_size[3];
  /** inverse cell size = \see DomainDecomposition::cell_size ^ -1. */
  double inv_cell_size[3];
  /** Array containing information about the interactions between the cells. */
  IA_Neighbor_List *cell_inter;
}  DomainDecomposition;


/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Information about the domain decomposition. */
extern DomainDecomposition dd;

/** Maximal skin size. This is a global variable wwhich can be read
    out by the user via \ref setmd in order to optimize the cell
    grid */
extern double max_skin;

/** Maximal number of cells per node. In order to avoid memory
 *  problems due to the cell grid one has to specify the maximal
 *  number of \ref cells::cells . The corresponding callback function
 *  is \ref max_num_cells_callback. If the number of cells \ref
 *  n_cells, is larger than max_num_cells the cell grid is
 *  reduced. max_num_cells has to be larger than 27, e.g one inner
 *  cell.  max_num_cells is initialized with the default value
 *  specified in \ref config.h: \ref CELLS_MAX_NUM_CELLS.
 */
extern int max_num_cells;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Re-derives the topology dimensions after the NpT-integrator
    changed the box-length. Note that no changes occur to the
    cell structure itself, use \ref dd_create_cell_grid for that.
    @param scal1 isotropic scaling factor by which each \ref box_l[i] changed. */
void dd_NpT_update_cell_grid(double scal1);

/** Initialize the topology. The argument is a list of cell pointers,
    containing particles that have to be sorted into new cells. The
    particles might not belong to this node.  This procedure is used
    when particle data or cell structure has changed and the cell
    structure has to be reinitialized. This also includes setting up
    the cell_structure array.
    @param cl List of cell pointers with particles to be stored in the
    new cell system.
*/
void dd_topology_init(CellPList *cl);

/** Called when the current cell structure is invalidated because for
    example the box length has changed. This procedure may NOT destroy
    the old inner and ghost cells, but it should free all other
    organizational data. Note that parameters like the box length or
    the node_grid may already have changed. Therefore organizational
    data has to be stored independently from variables that may be
    changed from outside. */
void dd_topology_release();

/** Just resort the particles. Used during integration. The particles
    are stored in the cell structure.

    @param global_flag Use DD_GLOBAL_EXCHANGE for global exchange and
    DD_NEIGHBOR_EXCHANGE for neighbor exchange (recommended for use within
    Molecular dynamics, or any other integration scheme using only local
    particle moves) 
*/
void dd_exchange_and_sort_particles(int global_flag);

/** implements \ref CellStructure::position_to_cell. */
Cell *dd_position_to_cell(double pos[3]);

/** Callback for setmd maxnumcells (maxnumcells >= 27). 
    see also \ref max_num_cells */
int max_num_cells_callback(Tcl_Interp *interp, void *_data);

/*@}*/

#endif
