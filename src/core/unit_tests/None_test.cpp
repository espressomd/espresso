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
 * Unit tests for the ScriptInterface::None class.
 *
 */

#define BOOST_TEST_MODULE None test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "None.hpp"
using ScriptInterface::None;

BOOST_AUTO_TEST_CASE(constructor_bool) {
  static_assert(!None{}, "");
  static_assert(!None{nullptr}, "");
}

BOOST_AUTO_TEST_CASE(comparison) {
  static_assert(None{} == None{}, "");
  static_assert(!(None{} != None{}), "");
  static_assert(!(None{} < None{}), "");
}

BOOST_AUTO_TEST_CASE(from_nullptr) {
  auto return_nullptr = []() -> None { return nullptr; };
  BOOST_CHECK(!return_nullptr());
}
