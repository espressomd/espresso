/*
  Copyright (C) 2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::Span test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/DeviceArray.hpp"
using Utils::DeviceArray;

#include <numeric>
#include <vector>
#include <array>

BOOST_AUTO_TEST_CASE(const_expr_ctor) {
  static_assert(4 == DeviceArray<int, 4>().size(), "");
  static_assert(4 == DeviceArray<int, 4>().max_size(), "");
}

BOOST_AUTO_TEST_CASE(array_ctor) {
  DeviceArray<int, 4> a;
  DeviceArray<int, 0> b;

  BOOST_CHECK_EQUAL(a.size(), 4);
  BOOST_CHECK_EQUAL(a.max_size(), 4);
  BOOST_CHECK_EQUAL(b.size(), 0);
  BOOST_CHECK_EQUAL(b.max_size(), 0);
}

BOOST_AUTO_TEST_CASE(iterators) {
  auto a = DeviceArray<int, 4>{{1, 2, 3, 4}};

  BOOST_CHECK(*(a.begin()) == 1);
  BOOST_CHECK(*(a.cbegin()) == 1);
  BOOST_CHECK(*(a.end()- 1) == 4);
  BOOST_CHECK(*(a.cend() -1) == 4);
}

BOOST_AUTO_TEST_CASE(element_access) {
  auto a = DeviceArray<int, 5>{{5,6,7,8,9}};

  int c = 5;
  for (DeviceArray<int, 5>::size_type i = 0; i < a.size(); ++i) {
    BOOST_CHECK(a.at(i) == c);
    BOOST_CHECK(a[i] == c);
    ++c;
  }

  BOOST_CHECK_THROW(a.at(a.size()), std::out_of_range);
}
