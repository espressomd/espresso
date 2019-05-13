/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef UTILS_COORDINATE_TRANSFORMATION_HPP
#define UTILS_COORDINATE_TRANSFORMATION_HPP

#include "utils/Vector.hpp"
#include "utils/constants.hpp"
#include "vec_rotate.hpp"

namespace Utils {

/** \brief Transform the given 3D position to cylinder coordinates with
 * longitudinal axis aligned with axis parameter.
 */
inline Vector3d
transform_coordinate_cartesian_to_cylinder(const Vector3d &pos,
                                           const std::string &axis) {
  static const Vector3d x_axis{1.0, 0.0, 0.0};
  static const Vector3d y_axis{0.0, 1.0, 0.0};
  Vector3d rotated_pos = pos;
  if (axis == "x") {
    rotated_pos = vec_rotate(y_axis, -Utils::pi() / 2.0, pos);
  } else if (axis == "y") {
    rotated_pos = vec_rotate(x_axis, Utils::pi() / 2.0, pos);
  }
  double r = std::sqrt(rotated_pos[0] * rotated_pos[0] +
                       rotated_pos[1] * rotated_pos[1]);
  double phi = std::atan2(rotated_pos[1], rotated_pos[0]);
  return Vector3d{r, phi, rotated_pos[2]};
}

/** \brief Transform the given 3D vector to cylinder coordinates with
 * longitudinal axis aligned with axis parameter.
 */
inline Vector3d transform_vector_cartesian_to_cylinder(
    Vector3d const &vec, std::string const &axis, Vector3d const &pos) {
  double v_r, v_phi, v_z;
  static const Vector3d x_axis{1.0, 0.0, 0.0};
  static const Vector3d y_axis{0.0, 1.0, 0.0};
  Vector3d rotated_vec = vec;
  Vector3d rotated_pos = pos;
  if (axis == "x") {
    rotated_vec = vec_rotate(y_axis, -Utils::pi() / 2.0, vec);
    rotated_pos = vec_rotate(y_axis, -Utils::pi() / 2.0, pos);
  } else if (axis == "y") {
    rotated_vec = vec_rotate(x_axis, Utils::pi() / 2.0, vec);
    rotated_pos = vec_rotate(x_axis, Utils::pi() / 2.0, pos);
  }
  // Coordinate transform the velocities.
  // v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
  v_r = (rotated_pos[0] * rotated_vec[0] + rotated_pos[1] * rotated_vec[1]) /
        std::sqrt(rotated_pos[0] * rotated_pos[0] +
                  rotated_pos[1] * rotated_pos[1]);
  // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
  v_phi = (rotated_pos[0] * rotated_vec[1] - rotated_pos[1] * rotated_vec[0]) /
          (rotated_pos[0] * rotated_pos[0] + rotated_pos[1] * rotated_pos[1]);
  // v_z = v_z
  v_z = rotated_vec[2];
  return Vector3d{v_r, v_phi, v_z};
}

/**
 * @brief Coordinate transformation from cylinder to cartesian coordinates.
 */
inline Vector3d
transform_coordinate_cylinder_to_cartesian(Vector3d const &p,
                                           const std::string &axis) {
  Vector3d transformed{{p[0] * std::cos(p[1]), p[0] * std::sin(p[1]), p[2]}};
  static const Vector3d x_axis{1.0, 0.0, 0.0};
  static const Vector3d y_axis{0.0, 1.0, 0.0};
  if (axis == "x") {
    transformed = vec_rotate(y_axis, Utils::pi() / 2.0, transformed);
  } else if (axis == "y") {
    transformed = vec_rotate(x_axis, -Utils::pi() / 2.0, transformed);
  }
  return transformed;
}

} // namespace Utils
#endif
