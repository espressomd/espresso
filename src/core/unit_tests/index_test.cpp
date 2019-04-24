/*
  Copyright (C) 2017-2018 The ESPResSo project

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

/** \file int_pow_test.cpp Unit tests for the
 * Utils::ravel_index and Utils::unravel_index functions.
 */

#define BOOST_TEST_MODULE Utils::ravel_index test
#define BOOST_TEST_DYN_LINK
#include <array>
#include <boost/test/unit_test.hpp>

#include "utils/index.hpp"

BOOST_AUTO_TEST_CASE(ravel_index_test) {
  const std::array<std::size_t, 4> unravelled_indices{{12, 23, 5, 51}};
  const std::array<std::size_t, 4> dimensions{{15, 35, 6, 52}};
  auto const result = Utils::ravel_index(unravelled_indices, dimensions);
  BOOST_CHECK(result == 138527);
}

BOOST_AUTO_TEST_CASE(unravel_index_test) {
  const std::array<std::size_t, 4> dimensions{{15, 35, 6, 52}};
  const std::size_t ravelled_index = 138527;
  std::array<std::size_t, 4> result{};
  Utils::unravel_index(dimensions.begin(), dimensions.end(), result.begin(),
                       ravelled_index);
  BOOST_CHECK((result == std::array<std::size_t, 4>{{12, 23, 5, 51}}));
}

BOOST_AUTO_TEST_CASE(get_linear_index) {
  using Utils::get_linear_index;

  auto const grid_size = Utils::Vector3i{7, 5, 3};
  auto const index = Utils::Vector3i{5, 4, 2};

  /* 3 ints version */
  {
    auto const linear_index =
        get_linear_index(index[0], index[1], index[2], grid_size);
    BOOST_CHECK_EQUAL(linear_index, (5 + 4 * 7 + 2 * (7 * 5)));
  }

  /* Vector version */
  {
    BOOST_CHECK_EQUAL(get_linear_index(index[0], index[1], index[2], grid_size),
                      get_linear_index(index, grid_size));
  }
}