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
#ifndef _CELLS_H
#define _CELLS_H
/** \file
 *  This file contains everything related to the cell structure / cell
 *  system.
 *
 *  The cell system (\ref CellStructure) describes how particles are
 *  distributed on the cells and how particles of different cells
 *  (regardless if they reside on the same or different nodes)
 *  interact with each other. The following cell systems are implemented:
 *
 *  - domain decomposition: The simulation box is divided spatially
 *    into cells (see \ref domain_decomposition.hpp). This is suitable for
 *    short range interactions.
 *  - nsquare: The particles are distributed equally on all nodes
 *    regardless their spatial position (see \ref nsquare.hpp). This is
 *    suitable for long range interactions that cannot be treated by a
 *    special method like P3M (see \ref p3m.hpp).
 */

#include "CellStructure.hpp"

#include <utility>
#include <vector>

/** \name Flags for exchange_and_sort_particles: whether to do a global
 *  exchange or assume that particles did not move much (faster, used
 *  during integration, where moving far is a catastrophe anyways).
 */
/*@{*/

enum {
  /** Flag for exchange_and_sort_particles : Do neighbor exchange. */
  CELL_NEIGHBOR_EXCHANGE = 0,
  /** Flag for exchange_and_sort_particles : Do global exchange. */
  CELL_GLOBAL_EXCHANGE = 1
};

/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** list of all cells. */
extern std::vector<Cell> cells;

/** Type of cell structure in use. */
extern CellStructure cell_structure;

/** If non-zero, cell systems should reset the position for checking
 *  the Verlet criterion. Moreover, the Verlet list has to be rebuilt.
 */
extern bool rebuild_verletlist;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Reinitialize the cell structures.
 *  @param new_cs The new topology to use afterwards. May be set to
 *                @ref CELL_STRUCTURE_CURRENT for not changing it.
 *  @param range  Desired interaction range
 */
void cells_re_init(int new_cs, double range);

/** Sort the particles into the cells and initialize the ghost particle
 *  structures.
 */
void cells_resort_particles(int global_flag);

/** This function is called whenever the cell system has to be
 *  reinitialized, e.g. if cutoffs have changed, or the skin, grid, ...
 *  It calculates the maximal interaction range, and as said reinitializes
 *  the cells structure if something significant has changed.
 *
 *  If the fast flag is set, the routine should try to save time.
 *  Currently this means that if the maximal range decreased, it does
 *  not reorganize the particles. This is used in the NpT algorithm to
 *  avoid frequent reorganization of particles.
 *
 *  @param fast If true, do not try to optimize the cell size.
 */
void cells_on_geometry_change(bool fast);

/** Update ghost information. If needed,
 *  the particles are also resorted.
 */
void cells_update_ghosts(unsigned data_parts);

/** Calculate and return the total number of particles on this node. */
int cells_get_n_particles();

/**
 * @brief Get pairs closer than @p distance from the cells.
 *
 * Pairs are sorted so that first.id < second.id
 */
std::vector<std::pair<int, int>> mpi_get_pairs(double distance);

/** Check if a particle resorting is required. */
void check_resort_particles();

/*@}*/

/**
 * @brief Find the cell in which a particle is stored.
 *
 * Uses position_to_cell on p.r.p. If this is not on the node's domain,
 * uses position at last Verlet list rebuild (p.l.p_old).
 *
 * @return pointer to the cell or nullptr if the particle is not on the node
 */
Cell *find_current_cell(const Particle &p);

#endif
