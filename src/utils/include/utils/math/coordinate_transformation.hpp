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

/**
 * @file
 * Convert coordinates from the Cartesian system to the cylindrical system.
 * The transformation functions are provided with three overloads:
 * - one function for the trivial Cartesian <-> cylindrical transformation
 * - one function to transform from/to a cylindrical system with custom axis
 *   (extra @p axis argument, keep in mind the angle phi is under-defined)
 * - one function to transform from/to an oriented cylindrical system with
 *   custom axis (extra @p orientation argument, the angle phi is well-defined)
 */

#include "utils/Vector.hpp"
#include "utils/constants.hpp"
#include "utils/math/interval.hpp"
#include "utils/math/vec_rotate.hpp"

#include <cmath>

namespace Utils {

/**
 * @brief Coordinate transformation from Cartesian to cylindrical coordinates.
 * The origins and z-axis of the coordinate systems co-incide.
 * The @f$ \phi = 0 @f$ direction corresponds to the x-axis in the
 * original coordinate system.
 * @param pos    %Vector to transform
 */
inline Vector3d
transform_coordinate_cartesian_to_cylinder(Vector3d const &pos) {
  auto const r = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  auto const phi = std::atan2(pos[1], pos[0]);
  return {r, phi, pos[2]};
}

/**
 * @brief Coordinate transformation from Cartesian to cylindrical coordinates
 * with change of basis. The origins of the coordinate systems co-incide.
 *
 * If the parameter @p axis is not equal to <tt>[0, 0, 1]</tt>, the value
 * of the angle @f$ \phi @f$ in cylindrical coordinates is under-defined.
 * To fully define it, it is necessary to provide an orientation vector
 * in Cartesian coordinates that will be used as the reference point
 * (i.e. such that @f$ \phi = 0 @f$).
 *
 * If you don't need performance, you can call the overloaded function
 * @ref transform_coordinate_cartesian_to_cylinder(Vector3d const &, <!--
 * --> Vector3d const &, Vector3d const &) directly. Otherwise, call this
 * function once on the orientation vector and cache its @f$ \phi @f$ angle;
 * this angle can be later passed as argument @p phi0 for a 50-80% speed-up.
 *
 * @param pos    %Vector to transform
 * @param axis   Longitudinal axis of the cylindrical coordinates
 * @param phi0   Angle offset to get @f$ \phi = 0 @f$ (obtained by calling
 *               this function on the orientation vector instead of @p pos)
 */
inline Vector3d transform_coordinate_cartesian_to_cylinder(Vector3d const &pos,
                                                           Vector3d const &axis,
                                                           double phi0 = 0.) {
  static auto const z_axis = Vector3d{{0, 0, 1}};
  auto const angle = angle_between(axis, z_axis);
  auto const rotation_axis = Utils::vector_product(axis, z_axis).normalize();
  auto const pos_rotated = vec_rotate(rotation_axis, angle, pos);
  auto pos_cylinder = transform_coordinate_cartesian_to_cylinder(pos_rotated);
  if (phi0 != 0.) {
    pos_cylinder[1] = interval(pos_cylinder[1] - phi0, -pi(), pi());
  }
  return pos_cylinder;
}

/**
 * @brief Coordinate transformation from Cartesian to cylindrical coordinates
 * with change of basis. The origins of the coordinate systems co-incide.
 * @param pos    %Vector to transform
 * @param axis   Longitudinal axis of the cylindrical coordinates
 * @param orientation   Reference point (in untransformed coordinates) for
 *                      which @f$ \phi = 0 @f$
 */
inline Vector3d transform_coordinate_cartesian_to_cylinder(
    Vector3d const &pos, Vector3d const &axis, Vector3d const &orientation) {
  auto const phi0 =
      transform_coordinate_cartesian_to_cylinder(orientation, axis)[1];
  return transform_coordinate_cartesian_to_cylinder(pos, axis, phi0);
}

/**
 * @brief Coordinate transformation from cylindrical to Cartesian coordinates.
 * The origins and z-axis of the coordinate systems co-incide.
 * The @f$ \phi = 0 @f$ direction corresponds to the x-axis in the
 * transformed coordinate system.
 * @param pos    %Vector to transform
 */
inline Vector3d
transform_coordinate_cylinder_to_cartesian(Vector3d const &pos) {
  auto const &rho = pos[0];
  auto const &phi = pos[1];
  auto const &z = pos[2];
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

/**
 * @brief Coordinate transformation from cylindrical to Cartesian coordinates
 * with change of basis. The origins of the coordinate systems co-incide.
 *
 * If the parameter @p axis is not equal to <tt>[0, 0, 1]</tt>, the value
 * of the angle @f$ \phi @f$ in cylindrical coordinates is under-defined.
 * To fully define it, it is necessary to provide an orientation vector
 * in Cartesian coordinates that will be used as the reference point
 * (i.e. such that @f$ \phi = 0 @f$).
 *
 * If you don't need performance, you can call the overloaded function
 * @ref transform_coordinate_cylinder_to_cartesian(Vector3d const &, <!--
 * --> Vector3d const &, Vector3d const &) directly. Otherwise, call function
 * @ref transform_coordinate_cartesian_to_cylinder(Vector3d const &, <!--
 * --> Vector3d const &, double) once on the orientation vector and cache
 * its @f$ \phi @f$ angle; this angle can be later passed as argument
 * @p phi0 for an 80-100% speed-up.
 *
 * @param pos    %Vector to transform
 * @param axis   Longitudinal axis of the cylindrical coordinates
 * @param phi0   Angle offset to get @f$ \phi = 0 @f$
 */
inline Vector3d transform_coordinate_cylinder_to_cartesian(Vector3d const &pos,
                                                           Vector3d const &axis,
                                                           double phi0 = 0.) {
  static auto const z_axis = Vector3d{{0, 0, 1}};
  auto const angle = angle_between(z_axis, axis);
  auto const rotation_axis = Utils::vector_product(z_axis, axis).normalize();
  auto const pos_cyl = Vector3d{{pos[0], pos[1] + phi0, pos[2]}};
  auto const pos_cart = transform_coordinate_cylinder_to_cartesian(pos_cyl);
  auto const pos_rot = vec_rotate(rotation_axis, angle, pos_cart);
  return pos_rot;
}

/**
 * @brief Coordinate transformation from cylindrical to Cartesian coordinates
 * with change of basis. The origins of the coordinate systems co-incide.
 * @param pos    %Vector to transform
 * @param axis   Longitudinal axis of the cylindrical coordinates
 * @param orientation   Reference point (in Cartesian coordinates) for
 *                      which @f$ \phi = 0 @f$
 */
inline Vector3d transform_coordinate_cylinder_to_cartesian(
    Vector3d const &pos, Vector3d const &axis, Vector3d const &orientation) {
  auto const phi0 =
      transform_coordinate_cartesian_to_cylinder(orientation, axis)[1];
  return transform_coordinate_cylinder_to_cartesian(pos, axis, phi0);
}

/**
 * @brief Vector transformation from Cartesian to cylindrical coordinates.
 * @param vec    %Vector to transform
 * @param axis   Longitudinal axis of the cylindrical coordinates
 * @param pos    Origin of the vector
 */
inline Vector3d transform_vector_cartesian_to_cylinder(Vector3d const &vec,
                                                       Vector3d const &axis,
                                                       Vector3d const &pos) {
  static auto const z_axis = Vector3d{{0, 0, 1}};
  auto const angle = angle_between(axis, z_axis);
  auto const rotation_axis = Utils::vector_product(axis, z_axis).normalize();
  auto const rotated_pos = vec_rotate(rotation_axis, angle, pos);
  auto const rotated_vec = vec_rotate(rotation_axis, angle, vec);
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
