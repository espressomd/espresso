#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#define BOOST_TEST_MODULE link_cell test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/Cell.hpp"
#include "core/algorithm/verlet_ia.hpp"

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
    for (auto &n : cells) {
      if (&c != &n)
        c.m_neighbors.push_back(&n);
    }

    c.part = new Particle[n_part_per_cell];
    c.n = c.max = n_part_per_cell;

    for (unsigned i = 0; i < n_part_per_cell; ++i) {
      c.part[i].p.identity = id++;
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

  /* Now check the verlet lists */
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

  for (auto &c : cells) {
    delete[] c.part;
  }
}
