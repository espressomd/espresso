/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <utility>
#include <vector>

BOOST_AUTO_TEST_CASE(link_cell) {
  const unsigned n_cells = 10;
  const auto n_part_per_cell = 10;
  const auto n_part = n_cells * n_part_per_cell;

  std::vector<Cell> cells(n_cells);

  auto id = 0;
  for (auto &c : cells) {
    std::vector<Cell *> neighbors;

    for (auto &n : cells) {
      if (&c != &n)
        neighbors.push_back(&n);
    }

    c.m_neighbors = Neighbors<Cell *>(neighbors, {});

    c.particles().resize(n_part_per_cell);

    for (auto &p : c.particles()) {
      p.id() = id++;
    }
  }

  std::vector<std::pair<int, int>> lc_pairs;
  lc_pairs.reserve((n_part * (n_part - 1)) / 2);

  Algorithm::link_cell(cells.begin(), cells.end(),
                       [&lc_pairs](Particle const &p1, Particle const &p2) {
                         if (p1.id() <= p2.id())
                           lc_pairs.emplace_back(p1.id(), p2.id());
                       });

  BOOST_CHECK(lc_pairs.size() == (n_part * (n_part - 1) / 2));

  auto it = lc_pairs.begin();
  for (int i = 0; i < n_part; i++)
    for (int j = i + 1; j < n_part; j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }
}
