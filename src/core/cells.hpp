/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

/**
 * @file
 * This file contains everything related to the global cell structure / cell
 * system.
 *
 * The cell system (@ref CellStructure) describes how particles are
 * distributed on the cells and how particles of different cells
 * (regardless if they reside on the same or different nodes)
 * interact with each other. The following cell systems are implemented:
 *
 * - regular decomposition: The simulation box is divided spatially
 *   into cells (see @ref RegularDecomposition.hpp). This is suitable for
 *   short-range interactions.
 * - N-square: The particles are distributed equally on all nodes
 *   regardless their spatial position (see @ref AtomDecomposition.hpp).
 *   This is suitable for long-range interactions that cannot be treated
 *   by a special method like P3M.
 * - hybrid decomposition: Initializes both regular decomposition
 *   and N-square at the same time and is given a set of particle types
 *   @c n_square_types (see @ref HybridDecomposition.hpp). By default,
 *   particles will be distributed using the regular decomposition.
 *   For particles of the types defined as @c n_square_types the N-square
 *   method is used. This is suitable for systems containing lots of small
 *   particles with short-range interactions mixed with a few large
 *   particles with long-range interactions. There, the large particles
 *   should be treated using N-square.
 */

#ifndef ESPRESSO_SRC_CORE_CELLS_HPP
#define ESPRESSO_SRC_CORE_CELLS_HPP

#include "cell_system/Cell.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <utility>
#include <vector>

/** Flags for particle exchange and resorting: whether to do a global
 *  exchange or assume that particles did not move much (faster, used
 *  during integration, where moving far is a catastrophe anyways).
 */
enum {
  /** Do neighbor exchange. */
  CELL_NEIGHBOR_EXCHANGE = 0,
  /** Do global exchange. */
  CELL_GLOBAL_EXCHANGE = 1
};

/** Type of cell structure in use. */
extern CellStructure cell_structure;

/** Initialize cell structure @ref HybridDecomposition
 *  @param n_square_types   Types of particles to place in the N-square cells.
 *  @param cutoff_regular   Cutoff for the regular decomposition.
 */
void set_hybrid_decomposition(std::set<int> n_square_types,
                              double cutoff_regular);

/** Reinitialize the cell structures.
 *  @param new_cs The new topology to use afterwards.
 */
void cells_re_init(CellStructureType new_cs);

/** Update ghost information. If needed,
 *  the particles are also resorted.
 */
void cells_update_ghosts(unsigned data_parts);

/**
 * @brief Get pairs closer than @p distance from the cells.
 *
 * Pairs are sorted so that first.id < second.id
 */
std::vector<std::pair<int, int>> get_pairs(double distance);

/**
 * @brief Get pairs closer than @p distance if both their types are in @p types
 *
 * Pairs are sorted so that first.id < second.id
 */
std::vector<std::pair<int, int>>
get_pairs_of_types(double distance, std::vector<int> const &types);

/** Check if a particle resorting is required. */
void check_resort_particles();

/**
 * @brief Get ids of particles that are within a certain distance
 * of another particle.
 */
boost::optional<std::vector<int>> get_short_range_neighbors(int pid,
                                                            double distance);

/**
 * @brief Find the cell in which a particle is stored.
 *
 * Uses position_to_cell on p.pos(). If this is not on the node's domain,
 * uses position at last Verlet list rebuild (p.p_old()).
 *
 * @return pointer to the cell or nullptr if the particle is not on the node
 */
Cell *find_current_cell(Particle const &p);

class PairInfo {
public:
  PairInfo() = default;
  PairInfo(int _id1, int _id2, Utils::Vector3d _pos1, Utils::Vector3d _pos2,
           Utils::Vector3d _vec21, int _node)
      : id1(_id1), id2(_id2), pos1(_pos1), pos2(_pos2), vec21(_vec21),
        node(_node) {}
  int id1;
  int id2;
  Utils::Vector3d pos1;
  Utils::Vector3d pos2;
  Utils::Vector3d vec21;
  int node;
};

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, PairInfo &p, unsigned int const /* version */) {
  ar &p.id1;
  ar &p.id2;
  ar &p.pos1;
  ar &p.pos2;
  ar &p.vec21;
  ar &p.node;
}
} // namespace serialization
} // namespace boost

/**
 * @brief Returns pairs of particle ids, positions and distance as seen by the
 * non-bonded loop.
 */
std::vector<PairInfo> non_bonded_loop_trace(int rank);

#endif
