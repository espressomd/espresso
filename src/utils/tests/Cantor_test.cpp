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

#define BOOST_TEST_MODULE Utils::Hash::Cantor test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <utility>

#include "utils/hash/Cantor.hpp"

BOOST_AUTO_TEST_CASE(cantor_pairing) {
  Utils::Hash::CantorPairing cantor;
  BOOST_CHECK(cantor(std::make_pair(1, 0)) == cantor(std::make_pair(0, 1)));
}

BOOST_AUTO_TEST_CASE(cantor_compare) {
  Utils::Hash::CantorCompare cantor;
  BOOST_CHECK(cantor(std::make_pair(1, 0), std::make_pair(0, 1)) == true);
}
