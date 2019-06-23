/*
  Copyright (C) 2017-2018 The ESPResSo project

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
 * Unit tests for the Utils::for_each_pair.
 *
 */

#define BOOST_TEST_MODULE for_each_pair test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/for_each_pair.hpp"
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
