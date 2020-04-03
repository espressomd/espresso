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
#include <functional>
#include <vector>

#define BOOST_TEST_MODULE link_cell test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Cell.hpp"
#include "algorithm/verlet_ia.hpp"

void check_pairs(int n_part, std::vector<std::pair<int, int>> const &pairs) {
  BOOST_CHECK(pairs.size() == (n_part * (n_part - 1) / 2));

  auto it = pairs.begin();
  for (int i = 0; i < n_part; i++)
    for (int j = i + 1; j < n_part; j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }
}

/* Dummy distance */
struct Distance {
  bool interact;
};

/* Dummy interaction criterion */
struct VerletCriterion {
  bool operator()(Particle const &p1, Particle const &p2,
                  Distance const &d) const {
    return d.interact;
  }
};

BOOST_AUTO_TEST_CASE(verlet_ia) {
  const unsigned n_cells = 100;
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
      p.p.identity = id++;
    }
  }

  std::vector<std::pair<int, int>> pairs;
  pairs.reserve((n_part * (n_part - 1)) / 2);
  std::vector<unsigned> id_counts(n_part, 0u);

  /* Build VL */
  Algorithm::verlet_ia(
      cells.begin(), cells.end(),
      [&id_counts](Particle const &p) { id_counts[p.p.identity]++; },
      [&pairs](Particle const &p1, Particle const &p2, Distance const &) {
        pairs.emplace_back(p1.p.identity, p2.p.identity);
      },
      [](Particle const &p1, Particle const &p2) {
        return Distance{p1.p.identity <= p2.p.identity};
      },
      VerletCriterion{}, /* rebuild */ true);

  /* Check that the particle kernel has been executed exactly once for every
   * particle. */
  BOOST_CHECK(std::all_of(id_counts.begin(), id_counts.end(),
                          [](int count) { return count == 1; }));

  check_pairs(n_part, pairs);

  /* Reset everything */
  pairs.clear();
  std::fill(id_counts.begin(), id_counts.end(), 0);

  /* Now check the Verlet lists */
  Algorithm::verlet_ia(
      cells.begin(), cells.end(),
      [&id_counts](Particle const &p) { id_counts[p.p.identity]++; },
      [&pairs](Particle const &p1, Particle const &p2, Distance const &) {
        pairs.emplace_back(p1.p.identity, p2.p.identity);
      },
      [](Particle const &p1, Particle const &p2) {
        return Distance{p1.p.identity <= p2.p.identity};
      },
      VerletCriterion{}, /* rebuild */ false);

  /* Check that the particle kernel has been executed exactly once for every
   * particle. */
  BOOST_CHECK(std::all_of(id_counts.begin(), id_counts.end(),
                          [](int count) { return count == 1; }));

  check_pairs(n_part, pairs);

  /* Reset everything */
  pairs.clear();
  std::fill(id_counts.begin(), id_counts.end(), 0);

  /* Rebuild again */
  Algorithm::verlet_ia(
      cells.begin(), cells.end(),
      [&id_counts](Particle const &p) { id_counts[p.p.identity]++; },
      [&pairs](Particle const &p1, Particle const &p2, Distance const &) {
        pairs.emplace_back(p1.p.identity, p2.p.identity);
      },
      [](Particle const &p1, Particle const &p2) {
        return Distance{p1.p.identity <= p2.p.identity};
      },
      VerletCriterion{}, /* rebuild */ true);

  /* Check that the particle kernel has been executed exactly once for every
   * particle. */
  BOOST_CHECK(std::all_of(id_counts.begin(), id_counts.end(),
                          [](int count) { return count == 1; }));

  check_pairs(n_part, pairs);
}
