/*
  Copyright (C) 2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::as_const test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/as_const.hpp"
using Utils::as_const;

#include <type_traits>

BOOST_AUTO_TEST_CASE(non_const) {
  using expected = const volatile int &;
  using actual = decltype(as_const(std::declval<volatile int &>()));

  static_assert(std::is_same<actual, expected>::value, "");
}

BOOST_AUTO_TEST_CASE(const_) {
  using expected = const volatile int &;
  using actual = decltype(as_const(std::declval<const volatile int &>()));

  static_assert(std::is_same<actual, expected>::value, "");
}

BOOST_AUTO_TEST_CASE(value) {
  int i = 5;

  BOOST_CHECK(i == as_const(i));
}
