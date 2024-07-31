/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE make_lin_space test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/math/make_lin_space.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

BOOST_AUTO_TEST_CASE(make_lin_space_test) {
  using Utils::make_lin_space;
  constexpr auto tol = 100. * std::numeric_limits<double>::epsilon();

  /* With endpoint */
  {
    auto const start = 1.;
    auto const stop = 2.;
    auto const num = 13u;

    auto const lin_space =
        make_lin_space(start, stop, num, /* endpoint */ true);
    BOOST_CHECK_EQUAL(lin_space.size(), num);

    std::vector<double> values;
    std::ranges::copy(lin_space, std::back_inserter(values));
    BOOST_CHECK_EQUAL(values.front(), start);
    BOOST_CHECK_CLOSE(values.back(), stop, tol);

    auto const dx = (stop - start) / static_cast<double>(num - 1u);
    for (std::size_t i = 0u; i < values.size(); i++) {
      BOOST_CHECK_CLOSE(values.at(i), start + static_cast<double>(i) * dx, tol);
    }
  }

  /* Without endpoint */
  {
    auto const start = 1.;
    auto const stop = 2.;
    auto const num = 13u;

    auto const lin_space =
        make_lin_space(start, stop, num, /* endpoint */ false);
    BOOST_CHECK_EQUAL(lin_space.size(), num);

    std::vector<double> values;
    std::ranges::copy(lin_space, std::back_inserter(values));
    BOOST_CHECK_EQUAL(values.front(), start);
    BOOST_CHECK_LT(values.back(), stop);

    auto const dx = (stop - start) / static_cast<double>(num);
    for (std::size_t i = 0u; i < values.size(); i++) {
      BOOST_CHECK_CLOSE(values.at(i), start + static_cast<double>(i) * dx, tol);
    }
  }
}
