/*
 * Copyright (C) 2020 The ESPResSo project
 *
 * This file is part of ESPResSo.
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

#define BOOST_TEST_MODULE Utils::linear_interpolation test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/linear_interpolation.hpp"

#include <vector>

BOOST_AUTO_TEST_CASE(interpolation) {
  using Utils::linear_interpolation;
  // tabulated values for f(x) = x^2
  auto const table = std::vector<double>{{1., 4., 9.}};
  constexpr auto tol = 1e-12;
  BOOST_CHECK_SMALL(linear_interpolation(table, 1., 1., 1.) - 1., tol);
  BOOST_CHECK_SMALL(linear_interpolation(table, 1., 1., 2.) - 4., tol);
  BOOST_CHECK_SMALL(linear_interpolation(table, 1., 1., 1.5) - 2.5, tol);
  BOOST_CHECK_SMALL(linear_interpolation(table, 1., 1., 3. - 1e-12) - 9.,
                    1e-10);
}
