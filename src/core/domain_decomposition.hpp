/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

#include "BoxGeometry.hpp"
#include "Cell.hpp"
#include "LocalBox.hpp"
#include "ghosts.hpp"

#include "DomainDecomposition.hpp"

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Information about the domain decomposition. */
extern DomainDecomposition dd;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize the topology. The argument is a list of cell pointers,
 *  containing particles that have to be sorted into new cells. The
 *  particles might not belong to this node. This procedure is used
 *  when particle data or cell structure has changed and the cell
 *  structure has to be reinitialized. This also includes setting up
 *  the cell_structure array.
 *
 *  @param comm MPI communicator to use for the cell system.
 *  @param range Desired interaction range
 */
void dd_topology_init(const boost::mpi::communicator &comm, double range,
                      const BoxGeometry &box_geo,
                      const LocalBox<double> &local_geo);

void dd_exchange_and_sort_particles(int global, ParticleList *pl,
                                    std::vector<Cell *> &modified_cells);

/*@}*/

#endif
