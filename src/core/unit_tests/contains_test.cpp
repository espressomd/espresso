/*
  Copyright (C) 2019 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::contains test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/contains.hpp"

BOOST_AUTO_TEST_CASE(contains_test) {
  using Utils::contains;

  auto const search_value = 31;
  int const true_rng[] = {9, 10, search_value, 11};
  int const false_rng[] = {1, 2, 3};

  /* iterator version */
  {
    BOOST_CHECK(
        contains(std::begin(true_rng), std::end(true_rng), search_value));
    BOOST_CHECK(
        not contains(std::begin(false_rng), std::end(false_rng), search_value));
  }

  /* range version */
  {
    BOOST_CHECK(contains(true_rng, search_value));
    BOOST_CHECK(not contains(false_rng, search_value));
  }
}
