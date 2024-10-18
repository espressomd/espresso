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

#ifndef ESPRESSO_SRC_CORE_CELL_SYSTEM_REGULAR_DECOMPOSITION_HPP
#define ESPRESSO_SRC_CORE_CELL_SYSTEM_REGULAR_DECOMPOSITION_HPP

#include "cell_system/ParticleDecomposition.hpp"

#include "cell_system/Cell.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"
#include "ghosts.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/optional.hpp>

#include <vector>

/**
 * @brief Regular decomposition cell system.
 *
 * The domain of a node is split into a 3D cell grid with dimension
 * @ref RegularDecomposition::cell_grid "cell_grid". Together with one ghost
 * cell layer on each side the overall dimension of the ghost cell grid is
 * @ref RegularDecomposition::ghost_cell_grid "ghost_cell_grid". The regular
 * decomposition enables the use of the linked cell algorithm
 * which is in turn used for setting up the Verlet list for the
 * system. You can see a 2D graphical representation of the linked
 * cell grid below.
 *
 * \image html  linked_cells.gif "Linked cells structure"
 *
 * 2D representation of a linked cell grid:
 * <tt>cell_grid = {4,4}, ghost_cell_grid = {6,6}</tt>
 *
 * Each cell has @f$ 3^D @f$ neighbor cells. Since we deal with pair forces,
 * it is sufficient to calculate only half of the interactions (Newton's law:
 * action = reaction). We have chosen the upper half e.g. all neighbor
 * cells with a higher linear index (for cell 14 they are marked in light
 * blue). Caution: This implementation needs double sided ghost
 * communication! For single sided ghost communication one would need
 * some ghost-ghost cell interaction as well, which we do not need!
 */
struct RegularDecomposition : public ParticleDecomposition {
  /** Grid dimensions per node. */
  Utils::Vector3i cell_grid = {};
  /** Cell size. */
  Utils::Vector3d cell_size = {};
  /** Offset in global grid */
  Utils::Vector3i cell_offset = {};
  /** linked cell grid with ghost frame. */
  Utils::Vector3i ghost_cell_grid = {};
  /** inverse @ref RegularDecomposition::cell_size "cell_size". */
  Utils::Vector3d inv_cell_size = {};

  boost::mpi::communicator m_comm;
  BoxGeometry const &m_box;
  LocalBox<double> m_local_box;
  std::vector<Cell> cells;
  std::vector<Cell *> m_local_cells;
  std::vector<Cell *> m_ghost_cells;
  GhostCommunicator m_exchange_ghosts_comm;
  GhostCommunicator m_collect_ghost_force_comm;

public:
  RegularSingleCoreDecomposition(double range, BoxGeometry const &box_geo,
                                 LocalBox<double> const &local_geo);

  GhostCommunicator const &exchange_ghosts_comm() const override { return {}; }
  GhostCommunicator const &collect_ghost_force_comm() const override {
    return {};
  }

  Utils::Span<Cell *> local_cells() override {
    return Utils::make_span(m_local_cells);
  }
  Utils::Span<Cell *> ghost_cells() override {
    return Utils::make_span(m_ghost_cells);
  }

  /* Getter needed for HybridDecomposition */
  std::vector<Cell *> get_local_cells() const { return m_local_cells; }
  std::vector<Cell *> get_ghost_cells() const { return {}; };

  Cell *particle_to_cell(Particle const &p) override {
    return position_to_cell(p.pos());
  }

  void resort(bool global, std::vector<ParticleChange> &diff) override;
  Utils::Vector3d max_cutoff() const override;
  Utils::Vector3d max_range() const override;

  boost::optional<BoxGeometry> minimum_image_distance() const override {
    return {m_box};
  }

  BoxGeometry const &box() const override { return m_box; };

private:
  /** Fill @c m_local_cells list and @c m_ghost_cells list for use with regular
   *  decomposition.
   */
  void mark_cells();

  /** Fill a communication cell pointer list. Fill the cell pointers of
   *  all cells which are inside a rectangular subgrid of the 3D cell
   *  grid starting from the
   *  lower left corner @p lc up to the high top corner @p hc. The cell
   *  pointer list @p part_lists must already be large enough.
   *  \param part_lists  List of cell pointers to store the result.
   *  \param lc          lower left corner of the subgrid.
   *  \param hc          high up corner of the subgrid.
   */
  void fill_comm_cell_lists(ParticleList **part_lists,
                            Utils::Vector3i const &lc,
                            Utils::Vector3i const &hc);

  int calc_processor_min_num_cells() const;

  Cell *position_to_cell(const Utils::Vector3d &pos);

  /**
   * @brief Move particles into the cell system if it belongs to this node.
   *
   * Moves all particles from src into the local cell
   * system if they do belong here. Otherwise the
   * particles are moved into rest.
   *
   * @param src Particles to move.
   * @param rest Output list for left-over particles.
   * @param modified_cells Local cells that were touched.
   */
  void move_if_local(ParticleList &src, ParticleList &rest,
                     std::vector<ParticleChange> &modified_cells);

  /**
  /**
   *  @brief Calculate cell grid dimensions, cell sizes and number of cells.
   *
   *  Calculates the cell grid, based on the local box size and the range.
   *  If the number of cells is larger than @c max_num_cells,
   *  it increases @c max_range until the number of cells is
   *  smaller or equal to @c max_num_cells. It sets:
   *  @c cell_grid,
   *  @c ghost_cell_grid,
   *  @c cell_size, and
   *  @c inv_cell_size.
   *
   *  @param range interaction range. All pairs closer
   *               than this distance are found.
   */
  void create_cell_grid(double range);

  /** Init cell interactions for cell system regular decomposition.
   *  Initializes the interacting neighbor cell list of a cell.
   *  This list of interacting neighbor cells is used by the Verlet
   *  algorithm.
   */
  void init_cell_interactions();

  /** Maximal number of cells per node. In order to avoid memory
   *  problems due to the cell grid, one has to specify the maximal
   *  number of cells. If the number of cells is larger
   *  than @c max_num_cells, the cell grid is reduced.
   *  @c max_num_cells has to be larger than 27, e.g. one inner cell.
   */
  static constexpr int max_num_cells = 32768;
};

#endif
