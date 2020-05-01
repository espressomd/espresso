/*
 * Copyright (C) 2019 The ESPResSo project
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

#define BOOST_TEST_MODULE quaternion test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <limits>

#include "utils/math/quaternion.hpp"

constexpr int w = -5;
constexpr int x = 2;
constexpr int y = 3;
constexpr int z = 4;
using Vector4i = Utils::VectorXi<4>;
Vector4i const zero_quat{{0, 0, 0, 0}};
Vector4i const identity_quat{{1, 0, 0, 0}};
Vector4i const scalar_quat{{w, 0, 0, 0}};
Vector4i const vector_quat{{0, x, y, z}};
Vector4i const full_quat{{w, x, y, z}};

BOOST_AUTO_TEST_CASE(multiply_quaternions) {
  using Utils::multiply_quaternions;
  /* identities */
  BOOST_CHECK(multiply_quaternions(identity_quat, scalar_quat) == scalar_quat);
  BOOST_CHECK(multiply_quaternions(identity_quat, scalar_quat) == scalar_quat);
  BOOST_CHECK(multiply_quaternions(identity_quat, full_quat) == full_quat);
  BOOST_CHECK(multiply_quaternions(scalar_quat, scalar_quat) ==
              w * scalar_quat);
  BOOST_CHECK(multiply_quaternions(scalar_quat, vector_quat) ==
              w * vector_quat);
  BOOST_CHECK(multiply_quaternions(zero_quat, full_quat) == zero_quat);
  BOOST_CHECK(multiply_quaternions(vector_quat, vector_quat) ==
              -vector_quat.norm2() * identity_quat);
  /* other */
  Vector4i const reference_quat{{-4, -20, -30, -40}};
  BOOST_CHECK(multiply_quaternions(full_quat, full_quat) == reference_quat);
}

BOOST_AUTO_TEST_CASE(convert_quaternion_to_director) {
  using Utils::convert_quaternion_to_director;
  /* identities */
  BOOST_CHECK(convert_quaternion_to_director(identity_quat) ==
              (Utils::Vector3i{{0, 0, 1}}));
  BOOST_CHECK(convert_quaternion_to_director(scalar_quat) ==
              (Utils::Vector3i{{0, 0, w * w}}));
  BOOST_CHECK(
      convert_quaternion_to_director(vector_quat) ==
      (Utils::Vector3i{{2 * x * z, 2 * y * z, -x * x - y * y + z * z}}));
  /* other */
  Utils::Vector3i const reference_director{{-14, 44, 28}};
  BOOST_CHECK(convert_quaternion_to_director(full_quat) == reference_director);
}

BOOST_AUTO_TEST_CASE(convert_director_to_quaternion) {
  using Utils::convert_director_to_quaternion;
  using Utils::Vector3d;
  using Utils::Vector4d;
  double const cos_pi_4 = std::sqrt(2.) / 2.;
  constexpr double eps = std::numeric_limits<double>::epsilon();
#define CHECK_QUAT(input, ref)                                                 \
  BOOST_CHECK_LE((convert_director_to_quaternion(input) - ref).norm2(), eps);
  /* identities */
  CHECK_QUAT((Vector3d{{0, 0, 0}}), (Vector4d{{1, 0, 0, 0}}));
  CHECK_QUAT((Vector3d{{0, 0, +1}}), (Vector4d{{1, 0, 0, 0}}));
  CHECK_QUAT((Vector3d{{0, 0, -1}}), (Vector4d{{0, -1, 0, 0}}));
  CHECK_QUAT((Vector3d{{+1, 0, 0}}), (Vector4d{{+1, -1, +1, -1}} / 2.));
  CHECK_QUAT((Vector3d{{-1, 0, 0}}), (Vector4d{{-1, +1, +1, -1}} / 2.));
  CHECK_QUAT((Vector3d{{0, +1, 0}}), (Vector4d{{+1, -1, 0, 0}} * cos_pi_4));
  CHECK_QUAT((Vector3d{{0, -1, 0}}), (Vector4d{{0, 0, +1, -1}} * cos_pi_4));
  /* self-consistency */
  using Utils::convert_quaternion_to_director;
  for (int i = -2; i < 3; ++i) {
    for (int j = -2; j < 3; ++j) {
      for (int k = -2; k < 3; ++k) {
        if (i == 0 and j == 0 and k == 0) {
          continue;
        }
        // here the case where j is zero is also tested
        Vector3d input{{static_cast<double>(i), static_cast<double>(j),
                        static_cast<double>(k)}};
        Vector3d output = convert_quaternion_to_director(
            convert_director_to_quaternion(input));
        input.normalize();
        output.normalize();
        BOOST_CHECK_LE((input - output).norm2(), eps);
      }
    }
  }
}
