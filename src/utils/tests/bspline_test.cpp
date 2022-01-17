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

#define BOOST_TEST_MODULE bspline test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/math/bspline.hpp>

#include <boost/mpl/list.hpp>

#include <array>
#include <stdexcept>

template <typename T, T... values>
using integer_list = boost::mpl::list<boost::mpl::integral_c<T, values>...>;
using test_bspline_orders = integer_list<int, 1, 2, 3, 4, 5, 6, 7>;

BOOST_AUTO_TEST_CASE_TEMPLATE(bspline_normalization, T, test_bspline_orders) {
  // check that B-splines are normalized
  constexpr auto order = T::value;
  constexpr auto tol = 1e-10;
  constexpr std::array<double, 5> x_values{-0.49999, 0.25, 0., 0.25, 0.49999};

  for (auto const x : x_values) {
    double sum = 0;
    for (int i = 0; i < order; ++i) {
      sum += Utils::bspline<order>(i, x);
    }
    BOOST_CHECK_SMALL(sum - 1., tol);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(bspline_symmetry, T, test_bspline_orders) {
  // check that B-splines are symmetric
  constexpr auto order = T::value;
  constexpr auto order_mid = (order % 2 == 0) ? order / 2 : (order + 1) / 2;
  constexpr auto tol = 1e-10;
  constexpr std::array<double, 3> x_values{-0.49999, 0.25, 0.1};

  for (int i = 0; i < order_mid; ++i) {
    for (auto const x : x_values) {
      auto const b_left = Utils::bspline<order>(i, x);
      auto const b_right = Utils::bspline<order>(order - i - 1, -x);
      BOOST_CHECK_SMALL(b_left - b_right, tol);
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(bspline_derivatives, T, test_bspline_orders) {
  // check that B-splines derivatives are correct
  constexpr auto order = T::value;
  constexpr auto tol = 1e-8;
  constexpr std::array<double, 5> x_values{-0.49999, 0.25, 0., 0.25, 0.49999};

  // approximate a derivative using the two-point central difference formula
  auto bspline_d_approx = [](int i, double x, int order) {
    using Utils::bspline;
    constexpr auto h = 1e-6;
    return (bspline(i, x + h / 2, order) - bspline(i, x - h / 2, order)) / h;
  };

  for (int i = 0; i < order; ++i) {
    for (auto const x : x_values) {
      auto const d_val = Utils::bspline_d<order>(i, x);
      auto const d_ref = bspline_d_approx(i, x, order);
      BOOST_CHECK_SMALL(d_val - d_ref, tol);
    }
  }
}

BOOST_AUTO_TEST_CASE(exceptions) {
  BOOST_CHECK_THROW((Utils::bspline<2, double>(-100, 0.)), std::runtime_error);
  BOOST_CHECK_THROW((Utils::bspline_d<2, float>(-1, 0.f)), std::runtime_error);
}
