/*
 * Copyright (C) 2018-2019 The ESPResSo project
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

/* Unit test for Utils Type Traits. */

#define BOOST_TEST_MODULE type traits tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/type_traits.hpp>

#include <climits>

BOOST_AUTO_TEST_CASE(size_in_bits) {
  static_assert(CHAR_BIT == Utils::size_in_bits<char>::value, "");
  static_assert(CHAR_BIT * sizeof(int) == Utils::size_in_bits<int>::value, "");
  BOOST_TEST_PASSPOINT();
}
