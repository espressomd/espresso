/*
  Copyright (C) 2016-2018 The ESPResSo project

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

#include <utility>

#define BOOST_TEST_MODULE Range test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Range.hpp"
using Utils::Range;

BOOST_AUTO_TEST_CASE(values) {
  int a[10];

  auto r = Range<int *>(a, a + 10);

  BOOST_CHECK(r.begin() == a);
  BOOST_CHECK(r.end() == a + 10);
}

BOOST_AUTO_TEST_CASE(std_begin_end) {
  int a[10];

  auto r = Range<int *>(a, a + 10);

  BOOST_CHECK(std::begin(r) == r.begin());
  BOOST_CHECK(std::end(r) == r.end());
}

BOOST_AUTO_TEST_CASE(make_range) {
  int a[10];

  BOOST_CHECK(Range<int *>(a, a + 10) ==
              Utils::make_range(std::begin(a), std::begin(a) + 10));
}

BOOST_AUTO_TEST_CASE(empty) {
  int a[10];

  auto r = Range<int *>(std::begin(a), std::begin(a));
  BOOST_CHECK(r.empty());
  auto s = Range<int *>(std::begin(a), std::begin(a) + 10);
  BOOST_CHECK(!s.empty());
}

BOOST_AUTO_TEST_CASE(size) {
  int a[10];

  auto r = Range<int *>(std::begin(a), std::begin(a) + 10);
  BOOST_CHECK(r.size() == 10);
}

BOOST_AUTO_TEST_CASE(comparison) {
  int a[10];

  BOOST_CHECK(!(Range<int *>(std::begin(a), std::begin(a)) ==
                Range<int *>(std::begin(a), std::begin(a) + 10)));
  BOOST_CHECK(Range<int *>(a, a + 10) ==
              Range<int *>(std::begin(a), std::begin(a) + 10));
}
