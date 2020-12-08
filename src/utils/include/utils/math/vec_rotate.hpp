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
#ifndef UTILS_VEC_ROTATE_HPP
#define UTILS_VEC_ROTATE_HPP

#include "utils/Vector.hpp"
#include "utils/math/sqr.hpp"

#include <cmath>
#include <tuple>

namespace Utils {
/**
 * @brief Rotate a vector around an axis.
 *
 * Rodrigues' rotation formula:
 * @f$ \vec{v}_{\mathrm{rot}} =
 *      (\cos{\alpha})\vec{v} +
 *      (\sin{\alpha})\vec{k}\times\vec{v} +
 *      (1 - \cos{\alpha})(\vec{k}\cdot\vec{v})\vec{k} @f$
 *
 * @param axis The axis to rotate about
 * @param alpha Angle to rotate
 * @param vector %Vector to act on
 * @return Rotated vector
 */
inline Vector3d vec_rotate(const Vector3d &axis, double alpha,
                           const Vector3d &vector) {
  if (std::abs(alpha) > std::numeric_limits<double>::epsilon()) {
    auto const sina = std::sin(alpha);
    auto const cosa = std::cos(alpha);
    auto const a = Vector3d(axis).normalize();
    auto const &v = vector;
    return cosa * v + sina * vector_product(a, v) + (1. - cosa) * (a * v) * a;
  }
  return vector;
}

/**
 * @brief Determine rotation angle and axis for rotating vec onto target_vec.
 * @param vec Vector to be rotated
 * @param target_vec Target vector
 * @return rotation angle and rotation axis
 */
inline std::tuple<double, Vector3d>
rotation_params(Vector3d const &vec, Vector3d const &target_vec) {
  if (vec.normalized() != target_vec.normalized()) {
    auto const theta =
        std::acos(vec * target_vec / (vec.norm() * target_vec.norm()));
    auto const rotation_axis =
        Utils::vector_product(vec, target_vec).normalize();
    return std::make_tuple(theta, rotation_axis);
  }
  return std::make_tuple(0.0, Vector3d{});
}

} // namespace Utils

#endif
