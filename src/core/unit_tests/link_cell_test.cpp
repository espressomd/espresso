/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE link_cell test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "algorithm/link_cell.hpp"

#include "Particle.hpp"
#include "cell_system/Cell.hpp"
#include "particle_store/ParticleStore.hpp"

#include <Kokkos_Core.hpp>

#include <utility>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

BOOST_AUTO_TEST_CASE(link_cell) {
  auto const n_cells = 10u;
  auto const n_part_per_cell = 10u;
  auto const n_part = n_cells * n_part_per_cell;

  std::vector<Cell> cells(n_cells);

  // Cells hold store ROW indices and hand out views. Build a store with
  // n_part local rows (ids 0..n_part-1 in cell-major order), wire each cell
  // to it, and give each cell its consecutive row block.
  ParticleStore store{};
  {
    std::vector<Particle> seeds(n_part);
    store.mark_dirty();
    store.begin_rebuild(n_part, 0u);
    int seed_row = 0;
    for (auto &p : seeds) {
      store.assign_row(p, seed_row++);
    }
    store.finish_rebuild();
    for (int r = 0; r < static_cast<int>(n_part); ++r) {
      store.id(r) = r;
    }
  }

  int row = 0;
  for (auto &c : cells) {
    std::vector<Cell *> neighbors;

    for (auto &n : cells) {
      if (&c != &n)
        neighbors.push_back(&n);
    }

    c.m_neighbors = Neighbors<Cell *>(neighbors, {});

    c.set_store(store);
    // A cell's committed rows are a contiguous (offset, count) range.
    // Each cell owns its consecutive block [row, row + n_part_per_cell).
    c.set_range(static_cast<std::size_t>(row), n_part_per_cell);
    row += static_cast<int>(n_part_per_cell);
  }

  std::vector<std::pair<int, int>> lc_pairs;
  lc_pairs.reserve((n_part * (n_part - 1u)) / 2u);

  Algorithm::link_cell(cells.begin(), cells.end(),
                       [&lc_pairs](Particle const &p1, Particle const &p2) {
                         if (p1.id() <= p2.id())
                           lc_pairs.emplace_back(p1.id(), p2.id());
                       });

  BOOST_CHECK(lc_pairs.size() == (n_part * (n_part - 1u) / 2u));

  auto it = lc_pairs.begin();
  for (auto i = 0; i < static_cast<int>(n_part); i++)
    for (auto j = i + 1; j < static_cast<int>(n_part); j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }
}
