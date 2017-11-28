#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#define BOOST_TEST_MODULE link_cell test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/Cell.hpp"
#include "core/algorithm/link_cell.hpp"

BOOST_AUTO_TEST_CASE(link_cell) {
  const unsigned n_cells = 100;
  const auto n_part_per_cell = 10;
  const auto n_part = n_cells * n_part_per_cell;

  std::vector<Cell> cells(n_cells);

  auto id = 0;
  for (auto &c : cells) {
    for (auto &n : cells) {
      if (&c != &n)
        c.m_neighbors.push_back(std::ref(n));
    }

    c.part = new Particle[n_part_per_cell];
    c.n = c.max = n_part_per_cell;

    for (unsigned i = 0; i < n_part_per_cell; ++i) {
      c.part[i].p.identity = id++;
    }
  }

  std::vector<std::pair<int, int>> lc_pairs;
  lc_pairs.reserve((n_part * (n_part - 1)) / 2);
  std::vector<unsigned> id_counts(n_part, 0u);

  Algorithm::link_cell(
      cells.begin(), cells.end(),
      [&id_counts](Particle const &p) { id_counts[p.p.identity]++; },
      [&lc_pairs](Particle const &p1, Particle const &p2,
                  std::pair<int, int> d) {
        /* Check that the "distance function" has been called with the corect
         * arguments */
        BOOST_CHECK((d.first == p1.p.identity) && (d.second == p2.p.identity));
        if (p1.p.identity <= p2.p.identity)
          lc_pairs.emplace_back(p1.p.identity, p2.p.identity);
      },
      [](Particle const &p1, Particle const &p2) {
        return std::make_pair(p1.p.identity, p2.p.identity);
      });

  /* Check that the particle kernel has been executed exactly once for every
   * particle. */
  BOOST_CHECK(std::all_of(id_counts.begin(), id_counts.end(),
                          [](int count) { return count == 1; }));

  BOOST_CHECK(lc_pairs.size() == (n_part * (n_part - 1) / 2));

  auto it = lc_pairs.begin();
  for (int i = 0; i < n_part; i++)
    for (int j = i + 1; j < n_part; j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }

  for (auto &c : cells) {
    delete[] c.part;
  }
}
