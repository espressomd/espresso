/*
 * Copyright (C) 2022 The ESPResSo project
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

#define BOOST_TEST_MODULE special functions test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "electrostatics/specfunc.hpp"

#include <cmath>
#include <limits>

auto constexpr eps = 8. * 100. * std::numeric_limits<double>::epsilon();

BOOST_AUTO_TEST_CASE(hurwitz_zeta_function) {
  constexpr auto max_bits = 54.0;
  // test cases where an exact, closed-form expression exists
  auto delta = 0.025;
  auto x = delta / 2.;
  while (x < 0.25) {
    auto order = max_bits / 2. + 0.01;
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(hzeta(order, x), std::pow(x, -order), eps);
    x += delta;
  }
  x = delta / 2.;
  while (x < 1.) {
    auto order = max_bits + 0.01;
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(hzeta(order, x), std::pow(x, -order), eps);
    x += delta;
  }
  x = delta / 2.;
  while (x < 1.) {
    auto order = max_bits / 2. + 0.01;
    auto ref = std::pow(x, -order);
    ref *= (1. + std::pow(x / (1. + x), order) + std::pow(x / (2. + x), order));
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(hzeta(order, x), ref, eps);
    x += delta;
  }
  // Hurwitz zeta is a generalization of the Riemann zeta
  for (auto order = 2; order < 128; ++order) {
    BOOST_TEST_INFO("with parameter order = " << order);
    BOOST_CHECK_CLOSE(hzeta(order, 1.), std::riemann_zeta(order), 10. * eps);
  }
}

BOOST_AUTO_TEST_CASE(bessel_series_high_precision) {
  auto delta = 0.02;
  auto x = delta;
  while (x < 20.) {
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(K0(x), std::cyl_bessel_k(0, x), 10. * eps);
    BOOST_CHECK_CLOSE(K1(x), std::cyl_bessel_k(1, x), 10. * eps);
    x += delta;
  }
  delta = 0.5;
  while (x < 64.) {
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(K0(x), std::cyl_bessel_k(0, x), 10. * eps);
    BOOST_CHECK_CLOSE(K1(x), std::cyl_bessel_k(1, x), 10. * eps);
    x += delta;
  }
}

BOOST_AUTO_TEST_CASE(bessel_series_low_precision) {
  auto const check = [](double x, double tol) {
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(LPK0(x), std::cyl_bessel_k(0, x), tol);
    BOOST_TEST_INFO("with parameter x = " << x);
    BOOST_CHECK_CLOSE(LPK1(x), std::cyl_bessel_k(1, x), tol);
    BOOST_CHECK_CLOSE(LPK01(x).first, LPK0(x), eps);
    BOOST_CHECK_CLOSE(LPK01(x).second, LPK1(x), eps);
  };
  auto delta = 0.02;
  auto x = delta;
  // important breakpoints: x=2, x=8, x=23, x=27 (different interpolation)
  while (x < 2.) {
    check(x, 20. * eps);
    x += delta;
  }
  while (x < 8.) {
    check(x, 1e5 * eps);
    x += delta;
  }
  delta = 0.04;
  while (x < 13.) {
    check(x, 1e8 * eps);
    x += delta;
  }
  while (x < 23.) {
    check(x, 2e10 * eps);
    x += delta;
  }
  delta = 0.08;
  while (x < 27.) {
    check(x, 0.6);
    x += delta;
  }
  delta = 0.64;
  while (x < 64.) {
    check(x, 20.);
    x += delta;
  }
}
