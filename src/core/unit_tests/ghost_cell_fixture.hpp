/*
 * Copyright (C) 2026 The ESPResSo project
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

#pragma once

/** @file
 *  Store-backed cells for the ghost-communication unit tests.
 *
 *  A @ref Cell owns a contiguous range of @ref ParticleStore rows rather than
 *  a particle container, so a test that wants "a cell holding N particles"
 *  needs a store to back it. This fixture builds one store whose rows are
 *  handed out to cells in declaration order: the first @c n_local cells get
 *  local rows, the rest get ghost rows, matching the layout that
 *  @c ensure_particle_store_synchronized produces.
 */

#include "Particle.hpp"
#include "cell_system/Cell.hpp"
#include "particle_store/ParticleStore.hpp"

#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

namespace GhostTest {

class CellFixture {
  ParticleStore m_store;
  std::vector<Particle> m_seeds;
  std::vector<Cell> m_cells;

public:
  /**
   * @param counts   Particle count per cell, in store-row order.
   * @param n_local  Number of leading cells that hold LOCAL rows; the
   *                 remaining cells hold ghost rows.
   */
  CellFixture(std::vector<std::size_t> const &counts, std::size_t n_local) {
    auto const rows_in = [&counts](std::size_t first, std::size_t last) {
      return std::accumulate(
          counts.begin() + static_cast<std::ptrdiff_t>(first),
          counts.begin() + static_cast<std::ptrdiff_t>(last), std::size_t{0});
    };
    auto const n_local_rows = rows_in(0u, n_local);
    auto const n_total_rows = rows_in(0u, counts.size());

    m_seeds.resize(n_total_rows);
    m_store.mark_dirty();
    m_store.begin_rebuild(n_local_rows, n_total_rows - n_local_rows);
    int row = 0;
    for (auto &p : m_seeds) {
      m_store.assign_row(p, row++);
    }
    m_store.finish_rebuild();

    // Cells are laid out back to back, exactly as the permutation rebuild
    // leaves them.
    m_cells.resize(counts.size());
    std::size_t offset = 0u;
    for (std::size_t i = 0u; i < counts.size(); ++i) {
      m_cells[i].set_store(m_store);
      m_cells[i].set_range(offset, counts[i]);
      offset += counts[i];
    }
  }

  Cell &cell(std::size_t i) { return m_cells[i]; }
  Cell const &cell(std::size_t i) const { return m_cells[i]; }
  Cell *ptr(std::size_t i) { return &m_cells[i]; }

  /** @brief First particle view of cell @p i (the tests use one-row cells). */
  Particle front(std::size_t i) {
    return m_store.make_view(static_cast<int>(m_cells[i].offset()));
  }

  ParticleStore &store() { return m_store; }
};

} // namespace GhostTest
