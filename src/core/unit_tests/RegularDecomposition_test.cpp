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

#define BOOST_TEST_MODULE RegularDecomposition test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/RegularDecomposition.hpp"
#include "communication.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <optional>
#include <utility>
#include <vector>

BOOST_TEST_GLOBAL_FIXTURE(EspressoCoreGlobalConfig);

/* Red/black balance of the cell pairing: on a fully periodic global cell grid
 * with at least three cells per direction, every cell owns (red) exactly half
 * of its 26 neighbor cells, also across the periodic boundaries, and every
 * pair of neighboring cells is owned by exactly one of its two cells, on one
 * or several MPI ranks. */
BOOST_AUTO_TEST_CASE(red_black_balance) {
  auto const &node_grid = ::communicator.node_grid;
  ::communicator.set_node_grid(node_grid);
  // cells of size 1 (= interaction range), ceil(4 / node_grid) cells per
  // rank and direction: a 4x4x4 global cell grid on 1, 2 and 4 ranks per
  // direction, 6 cells per direction on 3 ranks
  Utils::Vector3i local_cells_per_dir;
  Utils::Vector3d box_l;
  for (auto i : {0u, 1u, 2u}) {
    local_cells_per_dir[i] = (4 + node_grid[i] - 1) / node_grid[i];
    box_l[i] = node_grid[i] * local_cells_per_dir[i];
  }
  BoxGeometry box;
  box.set_length(box_l);
  for (auto i : {0u, 1u, 2u}) {
    box.set_periodic(i, true);
  }
  auto const local_box = LocalBox::make_regular_decomposition(
      box.length(), ::communicator.calc_node_index(), node_grid);
  RegularDecomposition dd{::communicator.comm, 1., box, local_box,
                          std::nullopt};

  Utils::Vector3i global_grid;
  for (auto i : {0u, 1u, 2u}) {
    global_grid[i] = dd.cell_grid[i] * node_grid[i];
    BOOST_REQUIRE_EQUAL(dd.cell_grid[i], local_cells_per_dir[i]);
    BOOST_REQUIRE_GE(global_grid[i], 3);
  }
  auto const n_cells_global = global_grid[0] * global_grid[1] * global_grid[2];

  // folded global linear index of a cell (local or ghost) of this rank
  auto const global_index = [&](Cell const *cell) {
    auto const idx = static_cast<int>(cell - dd.cells.data());
    BOOST_REQUIRE(idx >= 0 and idx < static_cast<int>(dd.cells.size()));
    Utils::Vector3i const ghost_pos{
        idx % dd.ghost_cell_grid[0],
        (idx / dd.ghost_cell_grid[0]) % dd.ghost_cell_grid[1],
        idx / (dd.ghost_cell_grid[0] * dd.ghost_cell_grid[1])};
    Utils::Vector3i pos;
    for (auto i : {0u, 1u, 2u}) {
      pos[i] = (ghost_pos[i] - 1 + dd.cell_offset[i] + global_grid[i]) %
               global_grid[i];
    }
    return Utils::get_linear_index(pos, global_grid);
  };

  // unordered pairs as (smaller, larger) global index, flattened
  std::vector<int> local_pairs;
  for (auto cell : dd.local_cells()) {
    BOOST_CHECK_EQUAL(cell->neighbors().red().size(), 13u);
    BOOST_CHECK_EQUAL(cell->neighbors().black().size(), 13u);
    auto const ind1 = global_index(cell);
    for (auto neighbor : cell->neighbors().red()) {
      auto const ind2 = global_index(neighbor);
      BOOST_CHECK_NE(ind1, ind2);
      local_pairs.emplace_back(std::min(ind1, ind2));
      local_pairs.emplace_back(std::max(ind1, ind2));
    }
  }

  Utils::Mpi::gather_buffer(local_pairs, ::communicator.comm);
  if (::communicator.comm.rank() == 0) {
    std::map<std::pair<int, int>, int> counts;
    for (std::size_t k = 0; k < local_pairs.size(); k += 2) {
      ++counts[{local_pairs[k], local_pairs[k + 1]}];
    }
    // 13 pairs per cell (26 neighbors, each pair counted once), every pair
    // occurs exactly once
    BOOST_CHECK_EQUAL(counts.size(),
                      static_cast<std::size_t>(13 * n_cells_global));
    BOOST_CHECK(std::ranges::all_of(
        counts, [](auto const &kv) { return kv.second == 1; }));
  }
}
