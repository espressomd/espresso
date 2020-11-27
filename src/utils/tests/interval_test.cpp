/*
 * Copyright (C) 2020 The ESPResSo project
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

#define BOOST_TEST_MODULE interval test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/constants.hpp>
#include <utils/math/interval.hpp>

#include <cmath>

BOOST_AUTO_TEST_CASE(interval_test) {
  constexpr auto eps = 1e-14;
  constexpr auto pi = Utils::pi<double>();
  BOOST_CHECK_EQUAL(Utils::interval(-pi, -pi, pi), -pi);
  BOOST_CHECK_EQUAL(Utils::interval(pi, -pi, pi), pi);
  for (int i = -24; i < 24; ++i) {
    double const theta = i * pi / 6;
    double const wrapped = Utils::interval(theta, -pi, pi);
    double const ref_val = std::atan2(std::sin(theta), std::cos(theta));
    BOOST_CHECK_SMALL(std::abs(std::fmod(wrapped - ref_val, 2 * pi)), eps);
  }
}
