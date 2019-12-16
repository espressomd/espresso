/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef ESPRESSO_ROTATION_MATRIX_HPP
#define ESPRESSO_ROTATION_MATRIX_HPP

#include "utils/Vector.hpp"

namespace Utils {

/**
 *  @brief Rotation matrix for quaternion.
 *  Taken from @cite goldstein80a (Chapter 4.6 eq. (4.47)).
 *
 *  @param q Quaternion
 *  @return Corresponding rotation matrix
 */
template <class T> Matrix<T, 3, 3> rotation_matrix(Vector<T, 4> const &q) {
  auto const q0q0 = q[0] * q[0];
  auto const q1q1 = q[1] * q[1];
  auto const q2q2 = q[2] * q[2];
  auto const q3q3 = q[3] * q[3];

  // clang-format off
  return {
          {q0q0 + q1q1 - q2q2 - q3q3, 2 * (q[1] * q[2] - q[0] * q[3]), 2 * (q[1] * q[3] + q[0] * q[2])},
          {2 * (q[1] * q[2] + q[0] * q[3]), q0q0 - q1q1 + q2q2 - q3q3, 2 * (q[2] * q[3] - q[0] * q[1])},
          {2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), q0q0 - q1q1 - q2q2 + q3q3}
         };
  // clang-format on
}
} // namespace Utils

#endif // ESPRESSO_ROTATION_MATRIX_HPP
