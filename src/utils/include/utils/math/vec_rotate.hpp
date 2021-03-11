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

#include <boost/qvm/quat_vec_operations.hpp>

#include "utils/Vector.hpp"
#include "utils/math/sqr.hpp"
#include "utils/quaternion.hpp"

#include <cmath>
#include <tuple>

namespace Utils {
/**
 * @brief Rotate a vector around an axis.
 *
 * @param axis The axis to rotate about
 * @param angle Angle to rotate
 * @param vector %Vector to act on
 * @return Rotated vector
 */
inline Vector3d vec_rotate(const Vector3d &axis, double angle,
                           const Vector3d &vector) {
  if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    Quaternion<double> q = boost::qvm::rot_quat(axis, angle);
    return q * vector;
  }
  return vector;
}

/**
 * @brief Determine the angle between two vectors.
 */
inline double angle_between(Vector3d const &v1, Vector3d const &v2) {
  return std::acos(v1 * v2 / std::sqrt(v1.norm2() * v2.norm2()));
}

} // namespace Utils

#endif
