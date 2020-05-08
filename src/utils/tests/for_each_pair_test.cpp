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

#include <numeric>
#include <vector>

void check_pairs(int n_values, std::vector<std::pair<int, int>> &pairs) {
  /* Check number of pairs */
  BOOST_CHECK(pairs.size() == ((n_values - 1) * n_values) / 2);

  /* Check that all were visited in order */
  auto it = pairs.begin();
  for (int i = 0; i < n_values; i++)
    for (int j = i + 1; j < n_values; j++) {
      BOOST_CHECK((it->first == i) && (it->second == j));
      ++it;
    }
}

BOOST_AUTO_TEST_CASE(basic_check) {
  auto const n_values = 141;
  std::vector<int> vec(n_values);

  std::iota(vec.begin(), vec.end(), 0);

  std::vector<std::pair<int, int>> pairs;

  auto pair_inserter = [&pairs](int a, int b) { pairs.emplace_back(a, b); };

  /* Collect pairs */
  for_each_pair(vec.begin(), vec.end(), pair_inserter);

  /* Check the result */
  check_pairs(n_values, pairs);
}

BOOST_AUTO_TEST_CASE(range) {
  auto const n_values = 141;
  std::vector<int> vec(n_values);
  std::iota(vec.begin(), vec.end(), 0);

  std::vector<std::pair<int, int>> pairs;

  auto pair_inserter = [&pairs](int a, int b) { pairs.emplace_back(a, b); };

  /* Collect pairs */
  for_each_pair(vec, pair_inserter);

  /* Check the result */
  check_pairs(n_values, pairs);
}

void check_cartesian_pairs(int n_values_1, int n_values_2,
                           std::vector<std::pair<int, int>> &pairs,
                           bool exclude_diagonal = false) {
  /* Check number of pairs */
  BOOST_CHECK(pairs.size() ==
              n_values_1 * n_values_2 -
                  (exclude_diagonal ? std::min(n_values_1, n_values_2) : 0));

  /* Check that all were visited in order */
  auto it = pairs.begin();
  for (int i = 0; i < n_values_1; ++i) {
    for (int j = 0; j < n_values_2; ++j) {
      if (!(exclude_diagonal and i == j)) {
        BOOST_CHECK((it->first == i) && (it->second == j));
        ++it;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(basic_check_cartesian) {
  auto const n_values1 = 141;
  auto const n_values2 = 143;
  std::vector<int> vec1(n_values1);
  std::vector<int> vec2(n_values2);

  std::iota(vec1.begin(), vec1.end(), 0);
  std::iota(vec2.begin(), vec2.end(), 0);

  std::vector<std::pair<int, int>> pairs;

  auto pair_inserter = [&pairs](int a, int b) { pairs.emplace_back(a, b); };

  /* Collect pairs */
  for_each_cartesian_pair(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                          pair_inserter);

  /* Check the result */
  check_cartesian_pairs(n_values1, n_values2, pairs);
}

BOOST_AUTO_TEST_CASE(range_cartesian) {
  auto const n_values1 = 141;
  auto const n_values2 = 143;
  std::vector<int> vec1(n_values1);
  std::vector<int> vec2(n_values2);

  std::iota(vec1.begin(), vec1.end(), 0);
  std::iota(vec2.begin(), vec2.end(), 0);

  std::vector<std::pair<int, int>> pairs;

  auto pair_inserter = [&pairs](int a, int b) { pairs.emplace_back(a, b); };

  /* Collect pairs */
  for_each_cartesian_pair(vec1, vec2, pair_inserter);

  /* Check the result */
  check_cartesian_pairs(n_values1, n_values2, pairs);
}

BOOST_AUTO_TEST_CASE(basic_check_cartesian_if) {
  auto const n_values1 = 141;
  auto const n_values2 = 143;
  std::vector<int> vec1(n_values1);
  std::vector<int> vec2(n_values2);

  std::iota(vec1.begin(), vec1.end(), 0);
  std::iota(vec2.begin(), vec2.end(), 0);

  std::vector<std::pair<int, int>> pairs;

  auto pair_inserter = [&pairs](int a, int b) { pairs.emplace_back(a, b); };
  auto binary_comparison = [](int a, int b) { return a != b; };

  /* Collect pairs */
  for_each_cartesian_pair_if(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                             pair_inserter, binary_comparison);

  /* Check the result */
  check_cartesian_pairs(n_values1, n_values2, pairs, true);
}

BOOST_AUTO_TEST_CASE(range_cartesian_if) {
  auto const n_values1 = 141;
  auto const n_values2 = 143;
  std::vector<int> vec1(n_values1);
  std::vector<int> vec2(n_values2);

  std::iota(vec1.begin(), vec1.end(), 0);
  std::iota(vec2.begin(), vec2.end(), 0);

  std::vector<std::pair<int, int>> pairs;

  auto pair_inserter = [&pairs](int a, int b) { pairs.emplace_back(a, b); };
  auto binary_comparison = [](int a, int b) { return a != b; };

  /* Collect pairs */
  for_each_cartesian_pair_if(vec1, vec2, pair_inserter, binary_comparison);

  /* Check the result */
  check_cartesian_pairs(n_values1, n_values2, pairs, true);
}
