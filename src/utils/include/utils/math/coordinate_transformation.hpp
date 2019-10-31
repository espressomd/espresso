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
#ifndef UTILS_COORDINATE_TRANSFORMATION_HPP
#define UTILS_COORDINATE_TRANSFORMATION_HPP

#include "utils/Vector.hpp"
#include "utils/constants.hpp"
#include "utils/math/vec_rotate.hpp"

namespace Utils {

/** \brief Transform the given 3D position to cylinder coordinates with
 * longitudinal axis aligned with axis parameter.
 */
inline Vector3d
transform_coordinate_cartesian_to_cylinder(const Vector3d &pos,
                                           const Vector3d &axis) {
  static auto const z_axis = Vector3d{{0, 0, 1}};
  double theta;
  Vector3d rotation_axis;
  std::tie(theta, rotation_axis) = rotation_params(axis, z_axis);
  auto const rotated_pos = vec_rotate(rotation_axis, theta, pos);
  auto const r = std::sqrt(rotated_pos[0] * rotated_pos[0] +
                           rotated_pos[1] * rotated_pos[1]);
  auto const phi = std::atan2(rotated_pos[1], rotated_pos[0]);
  return Vector3d{r, phi, rotated_pos[2]};
}

/**
 * @brief Coordinate transformation from cylinder to cartesian coordinates.
 */
inline Vector3d
transform_coordinate_cylinder_to_cartesian(Vector3d const &pos,
                                           Vector3d const &axis) {
  Vector3d const transformed{
      {pos[0] * std::cos(pos[1]), pos[0] * std::sin(pos[1]), pos[2]}};
  static auto const z_axis = Vector3d{{0, 0, 1}};
  double theta;
  Vector3d rotation_axis;
  std::tie(theta, rotation_axis) = rotation_params(z_axis, axis);
  auto const rotated_pos = vec_rotate(rotation_axis, theta, transformed);
  return rotated_pos;
}

/** \brief Transform the given 3D vector to cylinder coordinates with
 * symmetry axis aligned with axis parameter.
 */
inline Vector3d transform_vector_cartesian_to_cylinder(Vector3d const &vec,
                                                       Vector3d const &axis,
                                                       Vector3d const &pos) {
  static auto const z_axis = Vector3d{{0, 0, 1}};
  double theta;
  Vector3d rotation_axis;
  std::tie(theta, rotation_axis) = rotation_params(axis, z_axis);
  auto const rotated_pos = vec_rotate(rotation_axis, theta, pos);
  auto const rotated_vec = vec_rotate(rotation_axis, theta, vec);
  // v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
  auto const v_r =
      (rotated_pos[0] * rotated_vec[0] + rotated_pos[1] * rotated_vec[1]) /
      std::sqrt(rotated_pos[0] * rotated_pos[0] +
                rotated_pos[1] * rotated_pos[1]);
  // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
  auto const v_phi =
      (rotated_pos[0] * rotated_vec[1] - rotated_pos[1] * rotated_vec[0]) /
      (rotated_pos[0] * rotated_pos[0] + rotated_pos[1] * rotated_pos[1]);
  return Vector3d{v_r, v_phi, rotated_vec[2]};
}

} // namespace Utils
#endif
