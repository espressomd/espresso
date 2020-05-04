/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#define BOOST_TEST_MODULE ParticleIterator test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "mock/Particle.hpp"
using Testing::Particle;

#include "mock/Cell.hpp"
using Cell = Testing::Cell<Particle>;

#include "ParticleIterator.hpp"

std::vector<std::unique_ptr<Cell>> make_cells(std::size_t n) {
  std::vector<std::unique_ptr<Cell>> cells(n);

  for (auto &c : cells) {
    c = std::make_unique<Cell>();
  }

  return cells;
}

BOOST_AUTO_TEST_CASE(completeness) {
  using cells_t = std::vector<std::unique_ptr<Cell>>;
  using iterator = ParticleIterator<typename cells_t::iterator>;

  auto const n_part = 123456;
  auto cells = make_cells(1000);

  /* Fill the cells */
  for (int i = 0; i < n_part; i++) {
    cells[i % cells.size()]->part.emplace_back(i);
  }

  std::vector<int> counts(n_part, 0);

  auto begin = iterator(cells.begin(), cells.end());
  auto end = iterator(cells.end());

  /* Iterator over parts and count occurrence */
  for (; begin != end; ++begin) {
    counts[begin->identity()]++;
  }

  /* Every particle should be visited exactly once. */
  BOOST_CHECK(
      std::all_of(counts.begin(), counts.end(), [](int i) { return i == 1; }));
}

BOOST_AUTO_TEST_CASE(skip_empty) {
  using cells_t = std::vector<std::unique_ptr<Cell>>;
  using iterator = ParticleIterator<typename cells_t::iterator>;
  auto cells = make_cells(3);

  cells[0]->part.emplace_back(0);
  cells[2]->part.emplace_back(1);

  auto begin = iterator(cells.begin(), cells.end());

  BOOST_CHECK(begin->identity() == 0);
  ++begin;
  BOOST_CHECK(begin->identity() == 1);
}

BOOST_AUTO_TEST_CASE(order) {
  using cells_t = std::vector<std::unique_ptr<Cell>>;
  using iterator = ParticleIterator<typename cells_t::iterator>;
  auto const n_cells = 10;

  auto cells = make_cells(n_cells);

  /* Fill the cells */
  for (int i = 0; i < n_cells; i++) {
    cells[i % cells.size()]->part.emplace_back(i);
  }

  auto begin = iterator(cells.begin(), cells.end());
  auto end = iterator(cells.end());

  std::vector<Particle> id_diff(n_cells, Particle{0});
  std::adjacent_difference(begin, end, id_diff.begin(),
                           [](Particle const &a, Particle const &b) {
                             return Particle{a.identity() - b.identity()};
                           });

  BOOST_CHECK(std::all_of(id_diff.begin() + 1, id_diff.end(),
                          [](Particle const &p) { return p.identity() == 1; }));
}

BOOST_AUTO_TEST_CASE(distance_overload) {
  using cells_t = std::vector<std::unique_ptr<Cell>>;
  using iterator = ParticleIterator<typename cells_t::iterator>;
  auto const n_cells = 10;

  auto cells = make_cells(n_cells);

  /* Fill the cells */
  for (int i = 0; i < n_cells; i++) {
    cells[i % cells.size()]->part.emplace_back(i);
  }

  auto begin = iterator(cells.begin(), cells.end());
  auto end = iterator(cells.end());

  BOOST_CHECK(distance(begin, end) == std::distance(begin, end));
  BOOST_CHECK(distance(begin, begin) == 0);
  BOOST_CHECK(distance(end, end) == 0);
}

BOOST_AUTO_TEST_CASE(empty_range_) {
  using cells_t = std::vector<std::unique_ptr<Cell>>;
  using iterator = ParticleIterator<typename cells_t::iterator>;

  auto cells = make_cells(0);

  auto begin = iterator(cells.begin(), cells.end());
  auto end = iterator(cells.end());

  BOOST_CHECK_EQUAL(distance(begin, end), 0);
}