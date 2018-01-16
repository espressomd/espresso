#ifndef UTILS_COORDINATE_TRANSFORMATION_HPP
#define UTILS_COORDINATE_TRANSFORMATION_HPP

#include "utils.hpp"
#include "global.hpp"

extern double time_step;

namespace Utils {

/** \brief Transform the given 3D position to cylinder coordinates with longitudinal axis aligned with axis parameter.
 */
inline ::Vector<3, double>
    transform_pos_to_cylinder_coordinates(const ::Vector<3, double> &pos, const std::string &axis) {
  static const ::Vector<3, double> x_axis{1.0, 0.0, 0.0};
  static const ::Vector<3, double> y_axis{0.0, 1.0, 0.0};
  ::Vector<3, double> rotated_pos;
  if (axis == "x") {
    rotated_pos = vec_rotate(y_axis, -PI/2.0, pos);
  } else if (axis == "y") {
    rotated_pos = vec_rotate(x_axis, PI/2.0, pos);
  }
  double r = std::sqrt(rotated_pos[0] * rotated_pos[0] + rotated_pos[1] * rotated_pos[1]);
  double phi = std::atan2(rotated_pos[1], rotated_pos[0]);
  return ::Vector<3, double>{r, phi, rotated_pos[2]};
}

/** \brief Transform the given 3D velocity to cylinder coordinates with longitudinal axis aligned with axis parameter.
 */ 
inline ::Vector<3, double>
    transform_vel_to_cylinder_coordinates(::Vector<3, double> const &vel, const std::string &axis, ::Vector<3, double> const &pos) {
  double v_r, v_phi, v_z;
  static const ::Vector<3, double> x_axis{1.0, 0.0, 0.0};
  static const ::Vector<3, double> y_axis{0.0, 1.0, 0.0};
  ::Vector<3, double> rotated_vel;
  if (axis == "x") {
    rotated_vel = vec_rotate(y_axis, -PI/2.0, vel);
  } else if (axis == "y") {
    rotated_vel = vec_rotate(x_axis, PI/2.0, vel);
  }
  // Coordinate transform the velocities and divide core velocities by
  // time_step to get MD units. v_r = (x * v_x + y * v_y) / sqrt(x^2 +
  // y^2)
  v_r = (pos[0] * rotated_vel[0] / time_step +
         pos[1] * rotated_vel[1] / time_step) /
         std::sqrt(pos[0] * pos[0] +
                   pos[1] * pos[1]);
  // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
  v_phi = (pos[0] * rotated_vel[1] / time_step -
           pos[1] * rotated_vel[0] / time_step) /
          (pos[0] * pos[0] +
           pos[1] * pos[1]);
  // v_z = v_z
  v_z = rotated_vel[2] / time_step;
  return ::Vector<3, double>{v_r, v_phi, v_z};
}

} //namespace Utils
#endif
