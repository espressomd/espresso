/*
 * Copyright (C) 2018-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE rotation_matrix test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/math/rotation_matrix.hpp>
#include <utils/math/vec_rotate.hpp>

using namespace Utils;

BOOST_AUTO_TEST_CASE(transpose_test) {
  // clang-format off
  auto const A = Matrix<VectorXi<2>, 3, 3>{
      {{0,0}, {0, 1}, {0, 2}},
      {{1,0}, {1, 1}, {1, 2}},
      {{2,0}, {2, 1}, {2, 2}}
  };

  auto const expected = Matrix<VectorXi<2>, 3, 3>{
      {{0, 0}, {1, 0}, {2, 0}},
      {{0, 1}, {1, 1}, {2, 1}},
      {{0, 2}, {1, 2}, {2, 2}}
  };
  // clang-format on

  auto const result = transpose(A);
  BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(rotation_matrix_test) {
  /* Check rotation matrix by comparing with vec_rotate
   * for given axis and angle. */

  using std::cos;
  using std::sin;

  auto const axis = Vector3d{1., 2., 3.}.normalize();
  auto const angle = 0.7;
  auto const q = Vector4d{cos(angle / 2), sin(angle / 2) * axis[0],
                          sin(angle / 2) * axis[1], sin(angle / 2) * axis[2]};

  auto const v = Vector3d{3., 2., 1.};
  auto const expected = vec_rotate(axis, -angle, v);
  auto const result = rotation_matrix(q) * v;
  BOOST_CHECK_SMALL((expected - result).norm2(),
                    std::numeric_limits<double>::epsilon());
}