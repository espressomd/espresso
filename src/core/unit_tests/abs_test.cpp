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
#define BOOST_TEST_MODULE abs test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/math/abs.hpp>

#include <cmath>

BOOST_AUTO_TEST_CASE(abs_test) {
  using Utils::abs;

  static_assert(std::is_same<float, decltype(abs(1.f))>::value, "");
  static_assert(std::is_same<double, decltype(abs(1.))>::value, "");

  BOOST_CHECK_EQUAL(std::abs(3.1415), abs(3.1415));
  BOOST_CHECK_EQUAL(std::abs(-3.1415), abs(-3.1415));
  BOOST_CHECK_EQUAL(std::abs(3.1415f), abs(3.1415f));
  BOOST_CHECK_EQUAL(std::abs(-3.1415f), abs(-3.1415f));
}