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

namespace Utils {

/** \brief Transform the given 3D position to cylinder coordinates with
 * longitudinal axis aligned with axis parameter.
 */
inline ::Vector<3, double>
transform_pos_to_cylinder_coordinates(const ::Vector<3, double> &pos,
                                      const std::string &axis) {
  static const ::Vector<3, double> x_axis{1.0, 0.0, 0.0};
  static const ::Vector<3, double> y_axis{0.0, 1.0, 0.0};
  ::Vector<3, double> rotated_pos = pos;
  if (axis == "x") {
    rotated_pos = vec_rotate(y_axis, -PI / 2.0, pos);
  } else if (axis == "y") {
    rotated_pos = vec_rotate(x_axis, PI / 2.0, pos);
  }
  double r = std::sqrt(rotated_pos[0] * rotated_pos[0] +
                       rotated_pos[1] * rotated_pos[1]);
  double phi = std::atan2(rotated_pos[1], rotated_pos[0]);
  return ::Vector<3, double>{r, phi, rotated_pos[2]};
}

/** \brief Transform the given 3D velocity to cylinder coordinates with
 * longitudinal axis aligned with axis parameter.
 */
inline ::Vector<3, double>
transform_vel_to_cylinder_coordinates(const ::Vector<3, double> &vel,
                                      const std::string &axis,
                                      const ::Vector<3, double> &pos) {
  double v_r, v_phi, v_z;
  static const ::Vector<3, double> x_axis{1.0, 0.0, 0.0};
  static const ::Vector<3, double> y_axis{0.0, 1.0, 0.0};
  ::Vector<3, double> rotated_vel = vel;
  ::Vector<3, double> rotated_pos = pos;
  if (axis == "x") {
    rotated_vel = vec_rotate(y_axis, -PI / 2.0, vel);
    rotated_pos = vec_rotate(y_axis, -PI / 2.0, pos);
  } else if (axis == "y") {
    rotated_vel = vec_rotate(x_axis, PI / 2.0, vel);
    rotated_pos = vec_rotate(x_axis, PI / 2.0, pos);
  }
  // Coordinate transform the velocities.
  // v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
  v_r = (rotated_pos[0] * rotated_vel[0] + rotated_pos[1] * rotated_vel[1]) /
        std::sqrt(rotated_pos[0] * rotated_pos[0] +
                  rotated_pos[1] * rotated_pos[1]);
  // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
  v_phi = (rotated_pos[0] * rotated_vel[1] - rotated_pos[1] * rotated_vel[0]) /
          (rotated_pos[0] * rotated_pos[0] + rotated_pos[1] * rotated_pos[1]);
  // v_z = v_z
  v_z = rotated_vel[2];
  return ::Vector<3, double>{v_r, v_phi, v_z};
}

} // namespace Utils
#endif
