#ifndef UTILS_VEC_ROTATE_HPP
#define UTILS_VEC_ROTATE_HPP

#include <cmath>
#include <tuple>

#include "utils/Vector.hpp"
#include "utils/math/sqr.hpp"

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

/**
 * @brief Determine rotation angle and axis for rotating vec onto target_vec.
 * @param vec Vector to be rotated
 * @param target_vec Target vector
 * @return rotation angle and rotation axis
 */
inline std::tuple<double, Vector3d>
rotation_params(Vector3d const &vec, Vector3d const &target_vec) {
  auto const theta =
      std::acos(vec * target_vec) / (vec.norm() * target_vec.norm());
  auto const rotation_axis = Utils::vector_product(vec, target_vec).normalize();
  return std::make_tuple(theta, rotation_axis);
}

} // namespace Utils

#endif
