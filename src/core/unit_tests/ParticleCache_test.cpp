/*
  Copyright (C) 2017-2018 The ESPResSo project
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Unit tests for the MpiCallbacks class.
 *
 */

#include <random>
#include <vector>

#include "ParticleCache.hpp"

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ParticleCache test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "mock/Particle.hpp"
#include <utils/List.hpp>

using Communication::MpiCallbacks;
namespace mpi = boost::mpi;

class Particle : public Testing::Particle {
public:
  Particle() = default;
  Particle(int id) : Testing::Particle(id) {}

  IntList bl;

  IntList &bonds() { return bl; }
  IntList const &bonds() const { return bl; }

  Particle flat_copy() const { return Particle(m_id); }

  template <typename Archive> void serialize(Archive &ar, unsigned int) {
    ar &m_id;
    ar &bl.n;
  }
};

using Particles = std::vector<Particle>;

void check_merge(unsigned size, unsigned split) {
  boost::container::flat_set<int> u, v;
  using value_type = typename boost::container::flat_set<int>::value_type;

  std::mt19937 mersenne_engine(size);
  std::uniform_int_distribution<int> dist;

  for (int i = 0; i < size; i++) {
    auto val = value_type{dist(mersenne_engine)};
    if (i % split)
      u.insert(val);
    else
      v.insert(val);
  }

  auto merge = detail::Merge<boost::container::flat_set<int>, std::less<>>{};
  auto w = merge(u, v);

  BOOST_CHECK(w.size() == u.size() + v.size());

  for (auto const &e : u) {
    BOOST_CHECK(w.find(e) != w.end());
  }

  for (auto const &e : v) {
    BOOST_CHECK(w.find(e) != w.end());
  }
}

BOOST_AUTO_TEST_CASE(detail_merge_equal) { check_merge(2000, 2); }
BOOST_AUTO_TEST_CASE(detail_merge_not_equal) { check_merge(2000, 3); }
BOOST_AUTO_TEST_CASE(detail_merge_empty_left) { check_merge(2000, 1); }
BOOST_AUTO_TEST_CASE(detail_merge_empty_right) { check_merge(2000, 2000); }

BOOST_AUTO_TEST_CASE(update) {
  Particles local_parts;
  mpi::communicator world;
  MpiCallbacks cb(world);

  auto const rank = cb.comm().rank();
  auto const size = cb.comm().size();
  auto const n_part = 10000;

  local_parts.reserve(n_part);

  for (int i = 0; i < n_part; i++) {
    local_parts.emplace_back(rank * n_part + i);
  }

  auto get_parts = [&local_parts]() -> Particles const & {
    return local_parts;
  };
  ParticleCache<decltype(get_parts)> part_cfg(cb, get_parts);

  if (rank == 0) {
    BOOST_CHECK(part_cfg.size() == size * n_part);

    for (int i = 0; i < size * n_part; i++) {
      BOOST_CHECK(i == part_cfg[i].identity());
    }
  } else
    cb.loop();
}

BOOST_AUTO_TEST_CASE(update_with_bonds) {
  auto const bond_lengths = std::array<int, 6>{{1, 2, 4, 9, 21, 0}};

  Particles local_parts;
  mpi::communicator world;
  MpiCallbacks cb(world);

  auto const rank = cb.comm().rank();
  auto const size = cb.comm().size();
  auto const n_part = 1234;

  local_parts.reserve(n_part);

  for (int i = 0; i < n_part; i++) {
    auto const id = rank * n_part + (n_part - i - 1);
    local_parts.emplace_back(id);
    auto const bond_length = bond_lengths[id % bond_lengths.size()];
    auto &part = local_parts.back();
    part.bl.resize(bond_length);

    std::fill(part.bl.begin(), part.bl.end(), id);
  }

  auto get_parts = [&local_parts]() -> Particles const & {
    return local_parts;
  };
  ParticleCache<decltype(get_parts)> part_cfg(cb, get_parts);

  if (rank == 0) {
    part_cfg.update_bonds();

    for (int i = 0; i < size * n_part; i++) {
      /* Check that the length is set correctly */
      BOOST_CHECK(part_cfg[i].bl.size() ==
                  bond_lengths[part_cfg[i].identity() % bond_lengths.size()]);
      /* Check that the content was copied correctly. */
      BOOST_CHECK(std::all_of(part_cfg[i].bl.begin(), part_cfg[i].bl.end(),
                              [&i](int j) { return j == i; }));
    }
  } else {
    cb.loop();
  }
}

BOOST_AUTO_TEST_CASE(iterators) {
  Particles local_parts;
  mpi::communicator world;
  MpiCallbacks cb(world);

  auto const rank = cb.comm().rank();
  auto const size = cb.comm().size();
  auto const n_part = 1000;

  local_parts.reserve(n_part);

  for (int i = 0; i < n_part; i++) {
    local_parts.emplace_back(rank * n_part + i);
  }

  auto get_parts = [&local_parts]() -> Particles const & {
    return local_parts;
  };

  ParticleCache<decltype(get_parts)> part_cfg(cb, get_parts);

  if (rank == 0) {
    BOOST_CHECK(part_cfg.size() == size * n_part);

    std::vector<int> id_counts(size * n_part, 0);
    for (auto &p : part_cfg) {
      id_counts[p.identity()]++;
    }

    /* Every id should have been visited exactly once... */
    BOOST_CHECK(std::all_of(id_counts.begin(), id_counts.end(),
                            [](int count) { return count == 1; }));
    /* and in the correct order. */
    BOOST_CHECK(std::is_sorted(part_cfg.begin(), part_cfg.end(),
                               [](Particle const &a, Particle const &b) {
                                 return a.identity() < b.identity();
                               }));
  } else {
    cb.loop();
  }
}

int main(int argc, char **argv) {
  mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
