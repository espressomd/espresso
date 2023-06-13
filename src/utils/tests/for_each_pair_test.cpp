/*
 * Copyright (C) 2017-2022 The ESPResSo project
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

#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <utility>
#include <vector>

using PairContainer = std::vector<std::pair<std::size_t, size_t>>;

struct F {
  std::vector<std::size_t> vec1;
  std::vector<std::size_t> vec2;
  PairContainer pairs;
  F()
      : vec1(std::vector<std::size_t>(14)), vec2(std::vector<std::size_t>(15)),
        pairs() {
    std::iota(vec1.begin(), vec1.end(), 0);
    std::iota(vec2.begin(), vec2.end(), 0);
  }
  auto pair_inserter() {
    return [&pairs = pairs](std::size_t a, std::size_t b) {
      pairs.emplace_back(a, b);
    };
  }
};

void check_pairs(std::size_t n_values, PairContainer const &pairs) {
  /* Check number of pairs */
  BOOST_CHECK(pairs.size() == ((n_values - 1) * n_values) / 2);

  /* Check that all were visited in order */
  auto it = pairs.begin();
  for (std::size_t i = 0; i < n_values; i++)
    for (std::size_t j = i + 1; j < n_values; j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }
}

BOOST_FIXTURE_TEST_CASE(basic_check, F) {
  /* Collect pairs */
  for_each_pair(vec1.begin(), vec1.end(), pair_inserter());

  /* Check the result */
  check_pairs(vec1.size(), pairs);
}

BOOST_FIXTURE_TEST_CASE(range, F) {
  /* Collect pairs */
  for_each_pair(vec1, pair_inserter());

  /* Check the result */
  check_pairs(vec1.size(), pairs);
}

void check_cartesian_pairs(std::size_t n_values_1, std::size_t n_values_2,
                           PairContainer const &pairs,
                           bool exclude_diagonal = false) {
  /* Check number of pairs */
  BOOST_CHECK(pairs.size() ==
              n_values_1 * n_values_2 -
                  (exclude_diagonal ? std::min(n_values_1, n_values_2) : 0));

  /* Check that all were visited in order */
  auto it = pairs.begin();
  for (std::size_t i = 0; i < n_values_1; ++i) {
    for (std::size_t j = 0; j < n_values_2; ++j) {
      if (!exclude_diagonal or i != j) {
        BOOST_CHECK((it->first == i) && (it->second == j));
        ++it;
      }
    }
  }
}

BOOST_FIXTURE_TEST_CASE(basic_check_cartesian, F) {
  /* Collect pairs */
  for_each_cartesian_pair(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                          pair_inserter());

  /* Check the result */
  check_cartesian_pairs(vec1.size(), vec2.size(), pairs);
}

BOOST_FIXTURE_TEST_CASE(range_cartesian, F) {
  /* Collect pairs */
  for_each_cartesian_pair(vec1, vec2, pair_inserter());

  /* Check the result */
  check_cartesian_pairs(vec1.size(), vec2.size(), pairs);
}

BOOST_FIXTURE_TEST_CASE(basic_check_cartesian_if, F) {
  /* Collect pairs */
  for_each_cartesian_pair_if(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                             pair_inserter(), std::not_equal_to<std::size_t>{});

  /* Check the result */
  check_cartesian_pairs(vec1.size(), vec2.size(), pairs, true);
}

BOOST_FIXTURE_TEST_CASE(range_cartesian_if, F) {
  /* Collect pairs */
  for_each_cartesian_pair_if(vec1, vec2, pair_inserter(),
                             std::not_equal_to<std::size_t>{});

  /* Check the result */
  check_cartesian_pairs(vec1.size(), vec2.size(), pairs, true);
}
