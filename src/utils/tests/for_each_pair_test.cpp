/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

/* Unit tests for the Utils::for_each_pair. */

#define BOOST_TEST_MODULE for_each_pair test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/for_each_pair.hpp"
using Utils::for_each_cartesian_pair;
using Utils::for_each_cartesian_pair_if;
using Utils::for_each_pair;

#include <functional>
#include <numeric>
#include <vector>

using PairContainer = std::vector<std::pair<size_t, size_t>>;

struct MyFixture {
  std::vector<size_t> vec1;
  std::vector<size_t> vec2;
  PairContainer pairs;
  MyFixture()
      : vec1(std::vector<size_t>(141)), vec2(std::vector<size_t>(143)),
        pairs() {
    std::iota(vec1.begin(), vec1.end(), 0);
    std::iota(vec2.begin(), vec2.end(), 0);
  }
  auto pair_inserter() {
    return [&pairs = pairs](size_t a, size_t b) { pairs.emplace_back(a, b); };
  }
};

void check_pairs(size_t n_values, PairContainer const &pairs) {
  /* Check number of pairs */
  BOOST_CHECK(pairs.size() == ((n_values - 1) * n_values) / 2);

  /* Check that all were visited in order */
  auto it = pairs.begin();
  for (size_t i = 0; i < n_values; i++)
    for (size_t j = i + 1; j < n_values; j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }
}

BOOST_AUTO_TEST_CASE(basic_check) {
  MyFixture f;

  /* Collect pairs */
  for_each_pair(f.vec1.begin(), f.vec1.end(), f.pair_inserter());

  /* Check the result */
  check_pairs(f.vec1.size(), f.pairs);
}

BOOST_AUTO_TEST_CASE(range) {
  MyFixture f;

  /* Collect pairs */
  for_each_pair(f.vec1, f.pair_inserter());

  /* Check the result */
  check_pairs(f.vec1.size(), f.pairs);
}

void check_cartesian_pairs(size_t n_values_1, size_t n_values_2,
                           PairContainer const &pairs,
                           bool exclude_diagonal = false) {
  /* Check number of pairs */
  BOOST_CHECK(pairs.size() ==
              n_values_1 * n_values_2 -
                  (exclude_diagonal ? std::min(n_values_1, n_values_2) : 0));

  /* Check that all were visited in order */
  auto it = pairs.begin();
  for (size_t i = 0; i < n_values_1; ++i) {
    for (size_t j = 0; j < n_values_2; ++j) {
      if (!(exclude_diagonal and i == j)) {
        BOOST_CHECK((it->first == i) && (it->second == j));
        ++it;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(basic_check_cartesian) {
  MyFixture f;

  /* Collect pairs */
  for_each_cartesian_pair(f.vec1.begin(), f.vec1.end(), f.vec2.begin(),
                          f.vec2.end(), f.pair_inserter());

  /* Check the result */
  check_cartesian_pairs(f.vec1.size(), f.vec2.size(), f.pairs);
}

BOOST_AUTO_TEST_CASE(range_cartesian) {
  MyFixture f;

  /* Collect pairs */
  for_each_cartesian_pair(f.vec1, f.vec2, f.pair_inserter());

  /* Check the result */
  check_cartesian_pairs(f.vec1.size(), f.vec2.size(), f.pairs);
}

BOOST_AUTO_TEST_CASE(basic_check_cartesian_if) {
  MyFixture f;

  /* Collect pairs */
  for_each_cartesian_pair_if(f.vec1.begin(), f.vec1.end(), f.vec2.begin(),
                             f.vec2.end(), f.pair_inserter(),
                             std::not_equal_to<size_t>{});

  /* Check the result */
  check_cartesian_pairs(f.vec1.size(), f.vec2.size(), f.pairs, true);
}

BOOST_AUTO_TEST_CASE(range_cartesian_if) {
  MyFixture f;

  /* Collect pairs */
  for_each_cartesian_pair_if(f.vec1, f.vec2, f.pair_inserter(),
                             std::not_equal_to<size_t>{});

  /* Check the result */
  check_cartesian_pairs(f.vec1.size(), f.vec2.size(), f.pairs, true);
}
