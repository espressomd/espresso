#ifndef ESPRESSO_ROTATION_MATRIX_HPP
#define ESPRESSO_ROTATION_MATRIX_HPP

#include "utils/Vector.hpp"

namespace Utils {

/**
 * @brief Rotation matrix for quaternion
    Taken from "Goldstein - Classical Mechanics" (Chapter 4.6 Eq. 4.47).

    @param q Quaternion

    @return Corresponding rotation matrix
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

/**
 * @brief Transpose matrix
 * @param m Input matrix
 * @return Transposed matrix
 */
template <class T> Matrix<T, 3, 3> transpose(Matrix<T, 3, 3> const &m) {
  return {
      {m[0][0], m[1][0], m[2][0]},
      {m[0][1], m[1][1], m[2][1]},
      {m[0][2], m[1][2], m[2][2]},
  };
}
} // namespace Utils

#endif // ESPRESSO_ROTATION_MATRIX_HPP
