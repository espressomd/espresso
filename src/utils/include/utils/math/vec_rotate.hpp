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
  auto const sina = std::sin(alpha);
  auto const cosa = std::cos(alpha);
  auto const a = Vector3d(axis).normalize();
  auto const &v = vector;
  return cosa * v + sina * vector_product(a, v) + (1. - cosa) * (a * v) * a;
}

/**
 * @brief Determine the angle between two vectors.
 */
inline double angle_between(Vector3d const &v1, Vector3d const &v2) {
  return std::acos(v1 * v2 / (v1.norm() * v2.norm()));
}

} // namespace Utils

#endif
