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

#include "ParticleDecomposition.hpp"

#include <boost/mpi/communicator.hpp>
#include <boost/range/numeric.hpp>
#include <utils/index.hpp>
#include <utils/mpi/cart_comm.hpp>

/** Structure containing the information about the cell grid used for domain
 *  decomposition.
 */
struct DomainDecomposition {
  DomainDecomposition() = default;

  /** Offset in global grid */
  Utils::Vector3i cell_offset = {};
  /** linked cell grid in nodes spatial domain. */
  Utils::Vector3i cell_grid = {};
  /** linked cell grid with ghost frame. */
  Utils::Vector3i ghost_cell_grid = {};
  /** cell size. */
  Utils::Vector3d cell_size = {};
  /** inverse cell size = \see DomainDecomposition::cell_size ^ -1. */
  Utils::Vector3d inv_cell_size = {};
  bool fully_connected[3] = {false, false, false};

  boost::mpi::communicator comm;
  BoxGeometry box_geo;
  LocalBox<double> local_geo;
  std::vector<Cell> cells;
  std::vector<Cell *> m_local_cells;
  std::vector<Cell *> m_ghost_cells;
  GhostCommunicator m_exchange_ghosts_comm;
  GhostCommunicator m_collect_ghost_force_comm;

  /** Fill local_cells list and ghost_cells list for use with domain
   *  decomposition.  \ref cells::cells is assumed to be a 3d grid with size
   *  \ref DomainDecomposition::ghost_cell_grid.
   */
  void mark_cells() {
    int cnt_c = 0;

    m_local_cells.clear();
    m_ghost_cells.clear();

    for (int o = 0; o < ghost_cell_grid[2]; o++)
      for (int n = 0; n < ghost_cell_grid[1]; n++)
        for (int m = 0; m < ghost_cell_grid[0]; m++) {
          if ((m > 0 && m < ghost_cell_grid[0] - 1 && n > 0 &&
               n < ghost_cell_grid[1] - 1 && o > 0 &&
               o < ghost_cell_grid[2] - 1))
            m_local_cells.push_back(&cells.at(cnt_c++));
          else
            m_ghost_cells.push_back(&cells.at(cnt_c++));
        }
  }

  /** Fill a communication cell pointer list. Fill the cell pointers of
   *  all cells which are inside a rectangular subgrid of the 3D cell
   *  grid (\ref DomainDecomposition::ghost_cell_grid) starting from the
   *  lower left corner lc up to the high top corner hc. The cell
   *  pointer list part_lists must already be large enough.
   *  \param part_lists  List of cell pointers to store the result.
   *  \param lc          lower left corner of the subgrid.
   *  \param hc          high up corner of the subgrid.
   */
  void fill_comm_cell_lists(ParticleList **part_lists,
                            Utils::Vector3i const &lc,
                            Utils::Vector3i const &hc) {
    for (int o = lc[0]; o <= hc[0]; o++)
      for (int n = lc[1]; n <= hc[1]; n++)
        for (int m = lc[2]; m <= hc[2]; m++) {
          auto const i = Utils::get_linear_index(o, n, m, ghost_cell_grid);

          *part_lists++ = &(cells.at(i).particles());
        }
  }

  /** @brief Maximal interaction range supported with
   *         the current geometry and node grid.
   * @return Per-direction maximal range.
   */
  Utils::Vector3d max_range() const {
    auto dir_max_range = [this](int i) {
      if (fully_connected[i]) {
        return std::numeric_limits<double>::infinity();
      }

      return std::min(0.5 * box_geo.length()[i], local_geo.length()[i]);
    };

    return {dir_max_range(0), dir_max_range(1), dir_max_range(2)};
  }

  int calc_processor_min_num_cells() const {
    /* the minimal number of cells can be lower if there are at least two nodes
       serving a direction,
       since this also ensures that the cell size is at most half the box
       length. However, if there is only one processor for a direction, there
       have to be at least two cells for this direction. */
    return boost::accumulate(Utils::Mpi::cart_get<3>(comm).dims, 1,
                             [](int n_cells, int grid) {
                               return (grid == 1) ? 2 * n_cells : n_cells;
                             });
  }

public:
  /** Maximal number of cells per node. In order to avoid memory
   *  problems due to the cell grid one has to specify the maximal
   *  number of cells. If the number of cells is larger
   *  than max_num_cells the cell grid is reduced.
   *  max_num_cells has to be larger than 27, e.g. one inner cell.
   */
  static constexpr int max_num_cells = 32768;
};

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

/** adjust the domain decomposition to a change in the geometry.
 *  Tries to speed up things if possible. If function will fail,
 *  in which case it returns false, and in this case no changes
 *  to the cell system are made. If the change succeeds, all
 *  pointers into the cell system stay valid.
 *
 *  @param fast If true do not optimize the cell size but
 *              return asap.
 *  @param range Desired interaction range
 *  @return If the change was possible.
 */
bool dd_on_geometry_change(bool fast, double range, const BoxGeometry &box_geo,
                           const LocalBox<double> &local_geo);

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

/** Just resort the particles. Used during integration. The particles
 *  are stored in the cell structure.
 *
 *  @param global Use DD_GLOBAL_EXCHANGE for global exchange and
 *      DD_NEIGHBOR_EXCHANGE for neighbor exchange (recommended for use within
 *      Molecular dynamics, or any other integration scheme using only local
 *      particle moves)
 *  @param pl     List of particles
 */
void dd_exchange_and_sort_particles(int global, ParticleList *pl,
                                    std::vector<Cell *> &modified_cells);

/*@}*/

#endif
