/*
 * Copyright (C) 2017-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE Utils::sinc test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/constants.hpp"
#include "utils/math/sinc.hpp"

#include <cmath>
#include <cstdlib>

BOOST_AUTO_TEST_CASE(zero) { BOOST_CHECK_EQUAL(Utils::sinc(0.0), 1.0); }

BOOST_AUTO_TEST_CASE(approx) {
  auto x = 0.001;
  while (x <= 0.11) {
    auto const approx = Utils::sinc(x);
    auto const pi_x = Utils::pi() * x;
    auto const exact = std::sin(pi_x) / (pi_x);
    BOOST_CHECK_SMALL(approx - exact, 1e-13);
    x += 0.01;
  }
}
