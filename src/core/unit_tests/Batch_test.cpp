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

/** \file NumeratedContainer_test.cpp Unit tests for the
 * Utils::NumeratedContainer class.
 *
 */

#define BOOST_TEST_MODULE Utils::Batch test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/utils/Batch.hpp"

BOOST_AUTO_TEST_CASE(batch) {
  int counter = 0;
  auto a = [&counter](int i) {
    BOOST_CHECK(i == 42);
    BOOST_CHECK(counter++ == 0);
  };
  auto b = [&counter](int i) {
    BOOST_CHECK(i == 42);
    BOOST_CHECK(counter++ == 1);
  };
  auto c = [&counter](int i) {
    BOOST_CHECK(i == 42);
    BOOST_CHECK(counter++ == 2);
  };

  auto batch = Utils::make_batch(a, b, c);

  batch(42);
  BOOST_CHECK(counter == 3);
}
