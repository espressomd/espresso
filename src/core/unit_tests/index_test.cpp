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

BOOST_AUTO_TEST_CASE(ravel_index) {
  const std::array<std::size_t, 4> unravelled_indices{{12, 23, 5, 51}};
  const std::array<std::size_t, 4> dimensions{{15, 35, 6, 52}};
  auto const result = Utils::ravel_index(unravelled_indices, dimensions);
  BOOST_CHECK(result == 138527);
}

BOOST_AUTO_TEST_CASE(unravel_index) {
  const std::array<std::size_t, 4> dimensions{{15, 35, 6, 52}};
  const std::size_t ravelled_index = 138527;
  BOOST_CHECK((Utils::unravel_index(dimensions, ravelled_index) ==
               std::vector<std::size_t>{{12, 23, 5, 51}}));
}
