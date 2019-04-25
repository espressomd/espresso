#ifndef UTILS_VEC_ROTATE_HPP
#define UTILS_VEC_ROTATE_HPP

#include "utils/Vector.hpp"
#include "utils/math/sqr.hpp"

#include <cmath>

namespace Utils {
/**
 * @brief Rotate a vector around an axis.
 *
 * @param axis The axis to rotate about
 * @param alpha Angle to rotate
 * @param vector Vector to act on
 * @return Rotated vector
 */
inline Vector3d vec_rotate(const Vector3d &axis, double alpha,
                           const Vector3d &vector) {
  auto const sina = std::sin(alpha);
  auto const cosa = std::cos(alpha);
  auto const a = Vector3d(axis).normalize();

  return {(cosa + sqr(a[0]) * (1 - cosa)) * vector[0] +
              (a[0] * a[1] * (1 - cosa) - a[2] * sina) * vector[1] +
              (a[0] * a[2] * (1 - cosa) + a[1] * sina) * vector[2],
          (a[0] * a[1] * (1 - cosa) + a[2] * sina) * vector[0] +
              (cosa + sqr(a[1]) * (1 - cosa)) * vector[1] +
              (a[1] * a[2] * (1 - cosa) - a[0] * sina) * vector[2],
          (a[0] * a[2] * (1 - cosa) - a[1] * sina) * vector[0] +
              (a[1] * a[2] * (1 - cosa) + a[0] * sina) * vector[1] +
              (cosa + sqr(a[2]) * (1 - cosa)) * vector[2]};
}

} // namespace Utils

#endif
