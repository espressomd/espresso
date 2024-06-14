/*
 * Copyright (C) 2017-2024 The ESPResSo project
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

#define BOOST_TEST_MODULE "Special mathematical functions"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "p3m/math.hpp"

#include <utils/math/sqr.hpp>

#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <stdexcept>

/** @brief Compute the n-th term in the Taylor series of @c tan(x). */
static auto taylor_series_tangent(int n) {
  auto const two_power_2n = std::pow(2., 2 * n);
  auto const b2n = boost::math::bernoulli_b2n<double>(n);
  auto const f2n = boost::math::factorial<double>(2 * n);
  return std::pow(-1., n - 1) * two_power_2n * (two_power_2n - 1.) * b2n / f2n;
}

/** @brief Compute the n-th term in the Taylor series of @c sinc(x). */
static auto taylor_series_sinc(int n) {
  return std::pow(-1., n) / std::tgamma(2 * n + 2);
}

BOOST_AUTO_TEST_CASE(abs_test) {
  static_assert(std::is_same_v<float, decltype(math::abs(1.f))>);
  static_assert(std::is_same_v<double, decltype(math::abs(1.))>);

  BOOST_CHECK_EQUAL(math::abs(+3.1415), std::abs(+3.1415));
  BOOST_CHECK_EQUAL(math::abs(-3.1415), std::abs(-3.1415));
  BOOST_CHECK_EQUAL(math::abs(+3.1415f), std::abs(+3.1415f));
  BOOST_CHECK_EQUAL(math::abs(-3.1415f), std::abs(-3.1415f));
}

BOOST_AUTO_TEST_CASE(sinc_test) {
  // edge case
  BOOST_CHECK_EQUAL(math::sinc(0.0), 1.0);

  // check against standard math functions
  auto x = 0.001;
  while (x <= 0.11) {
    auto const approx = math::sinc(x);
    auto const pi_x = std::numbers::pi * x;
    auto const exact = std::sin(pi_x) / (pi_x);
    BOOST_CHECK_SMALL(approx - exact, 1e-13);
    x += 0.01;
  }

  // check Taylor expansion
  auto const series = std::views::iota(0, 5) | std::views::reverse |
                      std::views::transform(taylor_series_sinc);
  for (auto const x : {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}) {
    auto const pix2 = Utils::sqr(std::numbers::pi * x);
    auto const ref = std::accumulate(
        series.begin(), series.end(), 0.,
        [pix2](auto const &acc, auto const &c) { return acc * pix2 + c; });
    BOOST_CHECK_SMALL(math::sinc(x) - ref, 1e-13);
  }
}

BOOST_AUTO_TEST_CASE(analytic_cotangent_sum_test) {
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();
  auto const kernel = [](int n, double mesh_i, int cao) {
    return math::get_analytic_cotangent_sum_kernel(cao)(n, mesh_i);
  };

  // edge case: at theta = 0, aliasing sums are unity
  for (auto const cao : std::ranges::iota_view{1, 8}) {
    BOOST_CHECK_CLOSE(kernel(0, 0., cao), 1., tol);
  }
  // edge case: at theta = pi / 2, aliasing sums are equal to the cao-th term
  // in the tangent Taylor series
  for (auto const cao : std::ranges::iota_view{1, 8}) {
    BOOST_CHECK_CLOSE(kernel(1, 0.5, cao), taylor_series_tangent(cao), tol);
  }

  // check assertion
  for (auto const invalid_cao : {-1, 0, 8}) {
    BOOST_CHECK_THROW(kernel(1, 0., invalid_cao), std::logic_error);
  }
}
