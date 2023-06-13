/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE p3m test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "p3m/common.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

BOOST_AUTO_TEST_CASE(calc_meshift_false) {
  std::array<std::vector<int>, 3> const ref = {
      {std::vector<int>{0}, std::vector<int>{0, 1, -2, -1},
       std::vector<int>{0, 1, 2, 3, -3, -2, -1}}};

  auto const mesh = Utils::Vector3i{{1, 4, 7}};
  auto const val = detail::calc_meshift(mesh, false);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < ref[i].size(); ++j) {
      BOOST_CHECK_EQUAL(val[i][j], ref[i][j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(calc_meshift_true) {
  std::array<std::vector<int>, 3> const ref = {
      {std::vector<int>{0}, std::vector<int>{0, 1, 0, -1},
       std::vector<int>{0, 1, 2, 0, -3, -2, -1}}};

  auto const mesh = Utils::Vector3i{{1, 4, 7}};
  auto const val = detail::calc_meshift(mesh, true);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < ref[i].size(); ++j) {
      BOOST_CHECK_EQUAL(val[i][j], ref[i][j]);
    }
  }
}

#if defined(P3M) || defined(DP3M)
BOOST_AUTO_TEST_CASE(analytic_cotangent_sum) {
  auto constexpr kernel = p3m_analytic_cotangent_sum;
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();

  // check only trivial cases
  for (auto const cao : {1, 2, 3, 4, 5, 6, 7}) {
    BOOST_CHECK_CLOSE(kernel(0, 0., cao), 1., tol);
  }
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 1), 1., tol);
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 2), 1. / 3., tol);
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 3), 2. / 15., tol);
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 4), 17. / 315., tol);
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 5), 62. / 2835., tol);
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 6), 1382. / 155925., tol);
  BOOST_CHECK_CLOSE(kernel(1, 0.5, 7), 21844. / 6081075., tol);

  // check assertion
  for (auto const invalid_cao : {-1, 0, 8}) {
    BOOST_CHECK_THROW(kernel(1, 0., invalid_cao), std::logic_error);
  }
}
#endif // defined(P3M) || defined(DP3M)
