/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

/* Unit tests for the Utils::int_pow function. */

#define BOOST_TEST_MODULE int_pow test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/int_pow.hpp"
#include "utils/math/sqr.hpp"
using Utils::int_pow;
using Utils::sqr;

#include <limits>

auto const eps = std::numeric_limits<double>::epsilon();

/* Check that it can be used in constexpr context */
static_assert((int_pow<11>(2.), true), "");

/* Check that it can be used in constexpr context */
static_assert((sqr(2.), true), "");

BOOST_AUTO_TEST_CASE(even) {
  const double x = 3.14159;

  BOOST_CHECK(1 == int_pow<0>(x));
  BOOST_CHECK_CLOSE(x * x, int_pow<2>(x), 100. * eps);
  BOOST_CHECK_CLOSE((x * x) * (x * x), int_pow<4>(x), 100. * eps);
}

BOOST_AUTO_TEST_CASE(odd) {
  const double x = 3.14159;

  BOOST_CHECK(x == int_pow<1>(x));
  BOOST_CHECK_CLOSE((x * x) * x, int_pow<3>(x), 100. * eps);
  BOOST_CHECK_CLOSE((x * x) * (x * x) * x, int_pow<5>(x), 100. * eps);
}

BOOST_AUTO_TEST_CASE(square) {
  BOOST_CHECK_EQUAL(int_pow<2>(3.1415), sqr(3.1415));
}
