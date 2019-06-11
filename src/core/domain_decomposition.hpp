/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _DOMAIN_DECOMPOSITION_H
#define _DOMAIN_DECOMPOSITION_H

/** \file
 *
 *  This file contains everything related to the cell system: domain
 *  decomposition.
 *
 *  The simulation box is split into spatial domains for each node
 *  according to a Cartesian node grid (\ref node_grid).
 *
 *  The domain of a node is split into a 3D cell grid with dimension
 *  \ref DomainDecomposition::cell_grid. Together with one ghost cell
 *  layer on each side the overall dimension of the ghost cell grid is
 *  \ref DomainDecomposition::ghost_cell_grid. The domain
 *  decomposition enables one the use of the linked cell algorithm
 *  which is in turn used for setting up the Verlet list for the
 *  system. You can see a 2D graphical representation of the linked
 *  cell grid below.
 *
 *  \image html  linked_cells.gif "Linked cells structure"
 *
 *  2D representation of a linked cell grid: cell_grid =
 *  {4,4}, ghost_cell_grid = {6,6}
 *
 * Each cell has 3^D neighbor cells (For cell 14 they are
 * marked). Since we deal with pair forces, it is sufficient to
 * calculate only half of the interactions (Newtons law: action =
 * reaction). We have chosen the upper half e.g. all neighbor cells with
 * a higher linear index (For cell 14 they are marked in light
 * blue). Caution: This implementation needs double sided ghost
 * communication! For single sided ghost communication one would need
 * some ghost-ghost cell interaction as well, which we do not need!
 *
 *  For more information on cells, see \ref cells.hpp.
 *
 *  Implementation in domain_decomposition.cpp.
 */

#include "cells.hpp"

/** Structure containing the information about the cell grid used for domain
 *  decomposition.
 */
struct DomainDecomposition {
  DomainDecomposition()
      : cell_grid{0, 0, 0}, ghost_cell_grid{0, 0, 0}, cell_size{0, 0, 0},
        inv_cell_size{0, 0, 0} {}
  /** linked cell grid in nodes spatial domain. */
  int cell_grid[3];
  /** linked cell grid with ghost frame. */
  int ghost_cell_grid[3];
  /** cell size.
   *  Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range);
   * \endverbatim
   */
  double cell_size[3];
  /** inverse cell size = \see DomainDecomposition::cell_size ^ -1. */
  double inv_cell_size[3];
  bool fully_connected[3];
};

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Information about the domain decomposition. */
extern DomainDecomposition dd;

/** Maximal skin size. This is a global variable which can be read
 *  out by the user via the TCL command setmd in order to optimize the
 *  cell grid
 */
extern double max_skin;

/** Maximal number of cells per node. In order to avoid memory
 *  problems due to the cell grid one has to specify the maximal
 *  number of \ref cells::cells. If the number of cells is larger
 *  than max_num_cells the cell grid is reduced.
 *  max_num_cells has to be larger than 27, e.g. one inner cell.
 *  max_num_cells is initialized with the default value
 *  specified in \ref config.hpp as \ref CELLS_MAX_NUM_CELLS.
 */
extern int max_num_cells;

/** Minimal number of cells per node. This is mainly to avoid excessively large
 *  numbers of particles per cell, which will result in really large Verlet
 *  lists and eventually crash Espresso.
 */
extern int min_num_cells;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** adjust the domain decomposition to a change in the geometry.
 *  Tries to speed up things if possible.
 *
 *  @param flags  A combination of \ref CELL_FLAG_FAST and \ref
 *                CELL_FLAG_GRIDCHANGED, see documentation of \ref
 *                cells_on_geometry_change.
 *  @param grid   Number of nodes in each spatial dimension.
 */
void dd_on_geometry_change(int flags, const Utils::Vector3i &grid);

/** Initialize the topology. The argument is a list of cell pointers,
 *  containing particles that have to be sorted into new cells. The
 *  particles might not belong to this node.  This procedure is used
 *  when particle data or cell structure has changed and the cell
 *  structure has to be reinitialized. This also includes setting up
 *  the cell_structure array.
 *  @param cl    List of cell pointers with particles to be stored in the
 *               new cell system.
 *  @param grid  Number of nodes in each spatial dimension.
 */
void dd_topology_init(CellPList *cl, const Utils::Vector3i &grid);

/** Called when the current cell structure is invalidated because for
 *  example the box length has changed. This procedure may NOT destroy
 *  the old inner and ghost cells, but it should free all other
 *  organizational data. Note that parameters like the box length or
 *  the node_grid may already have changed. Therefore organizational
 *  data has to be stored independently from variables that may be
 *  changed from outside.
 */
void dd_topology_release();

/** Just resort the particles. Used during integration. The particles
 *  are stored in the cell structure.
 *
 *  @param global Use DD_GLOBAL_EXCHANGE for global exchange and
 *      DD_NEIGHBOR_EXCHANGE for neighbor exchange (recommended for use within
 *      Molecular dynamics, or any other integration scheme using only local
 *      particle moves)
 *  @param pl     List of particles
 *  @param grid   Number of nodes in each spatial dimension
 */
void dd_exchange_and_sort_particles(int global, ParticleList *pl,
                                    const Utils::Vector3i &grid);

/** calculate physical (processor) minimal number of cells */
int calc_processor_min_num_cells(const Utils::Vector3i &grid);

/** Fill a communication cell pointer list. Fill the cell pointers of
 *  all cells which are inside a rectangular subgrid of the 3D cell
 *  grid (\ref DomainDecomposition::ghost_cell_grid) starting from the
 *  lower left corner lc up to the high top corner hc. The cell
 *  pointer list part_lists must already be large enough.
 *  \param part_lists  List of cell pointers to store the result.
 *  \param lc          lower left corner of the subgrid.
 *  \param hc          high up corner of the subgrid.
 */
int dd_fill_comm_cell_lists(Cell **part_lists, int const lc[3],
                            int const hc[3]);

/** Of every two communication rounds, set the first receivers to prefetch and
 *  poststore
 */
void dd_assign_prefetches(GhostCommunicator *comm);

/*@}*/

#endif
