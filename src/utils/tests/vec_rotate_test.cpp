/*
 * Copyright (C) 2018-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE Utils::vec_rotate test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/constants.hpp>
#include <utils/math/vec_rotate.hpp>
using Utils::vec_rotate;

#include <cmath>
#include <limits>

BOOST_AUTO_TEST_CASE(rotation) {
  using std::cos;
  using std::sin;

  /* Axis */
  auto const k = Utils::Vector3d{1, 2, 3}.normalize();
  /* Angle */
  auto const t = 1.23;
  /* Original vector */
  auto const v = Utils::Vector3d{2, 3, 4};

  /* Rodrigues' formula from wikipedia */
  auto const expected =
      cos(t) * v + sin(t) * vector_product(k, v) + (1. - cos(t)) * (k * v) * k;

  auto const is = vec_rotate(k, t, v);
  auto const rel_diff = (expected - is).norm() / expected.norm();

  BOOST_CHECK(rel_diff < std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(rotation_params) {
  Utils::Vector3d v1 = {1.0, 0.0, 0.0};
  Utils::Vector3d v2 = {1.0, 1.0, 0.0};

  double angle;
  Utils::Vector3d rotation_axis;
  std::tie(angle, rotation_axis) = Utils::rotation_params(v1, v2);
  BOOST_CHECK_CLOSE(angle, Utils::pi() / 4.0, 1e-7);
  BOOST_CHECK_SMALL((rotation_axis * v1), 1e-7);
  BOOST_CHECK_SMALL((rotation_axis * v2), 1e-7);
}
