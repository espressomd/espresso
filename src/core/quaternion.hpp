/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
#ifndef QUATERNION_H
#define QUATERNION_H
/** \file
 *  Quaternion algebra.
 */

#include "config.hpp"

#ifdef ROTATION

#include <utils/Vector.hpp>
#include <utils/constants.hpp>

/** Multiply two quaternions */
template <typename T1, typename T2, typename T3>
void multiply_quaternions(const T1 &a, const T2 &b, T3 &result) {
  // Formula from http://www.j3d.org/matrix_faq/matrfaq_latest.html
  result[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  result[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
  result[2] = a[0] * b[2] + a[2] * b[0] + a[3] * b[1] - a[1] * b[3];
  result[3] = a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1];
}

inline Utils::Vector4d multiply_quaternions(const Utils::Vector4d &q,
                                            const Utils::Vector4d &p) {
  Utils::Vector4d res;
  multiply_quaternions(q, p, res);

  return res;
}

inline void convert_quat_to_director(const Utils::Vector4d &quat,
                                     Utils::Vector3d &director) {
  /* director */
  director[0] = 2 * (quat[1] * quat[3] + quat[0] * quat[2]);
  director[1] = 2 * (quat[2] * quat[3] - quat[0] * quat[1]);
  director[2] = (quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] +
                 quat[3] * quat[3]);
}

inline Utils::Vector3d convert_quat_to_director(const Utils::Vector4d &q) {
  Utils::Vector3d res;
  convert_quat_to_director(q, res);

  return res;
}

/** Convert director to quaternions */
inline int convert_director_to_quat(const Utils::Vector3d &d, Utils::Vector4d &quat) {
  double theta2, phi2;

  // Calculate magnitude of the given vector
  auto const dm = d.norm();

  // The vector needs to be != 0 to be converted into a quaternion
  if (dm < ROUND_ERROR_PREC) {
    return 1;
  }
  // Calculate angles
  auto const d_xy = sqrt(d[0] * d[0] + d[1] * d[1]);
  // If dipole points along z axis:
  if (d_xy == 0) {
    // We need to distinguish between (0,0,d_z) and (0,0,d_z)
    if (d[2] > 0)
      theta2 = 0;
    else
      theta2 = Utils::pi() / 2.;
    phi2 = 0;
  } else {
    // Here, we take care of all other directions
    // Here we suppose that theta2 = 0.5*theta and phi2 = 0.5*(phi -
    // Utils::pi()/2), where theta and phi - angles are in spherical coordinates
    theta2 = 0.5 * acos(d[2] / dm);
    if (d[1] < 0)
      phi2 = -0.5 * acos(d[0] / d_xy) - Utils::pi() * 0.25;
    else
      phi2 = 0.5 * acos(d[0] / d_xy) - Utils::pi() * 0.25;
  }

  // Calculate the quaternion from the angles
  auto const cos_theta2 = cos(theta2);
  auto const sin_theta2 = sin(theta2);
  auto const cos_phi2 = cos(phi2);
  auto const sin_phi2 = sin(phi2);
  quat[0] = cos_theta2 * cos_phi2;
  quat[1] = -sin_theta2 * cos_phi2;
  quat[2] = -sin_theta2 * sin_phi2;
  quat[3] = cos_theta2 * sin_phi2;

  return 0;
}

inline void normalize_quaternion(double *q) {
  double tmp = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  q[0] /= tmp;
  q[1] /= tmp;
  q[2] /= tmp;
  q[3] /= tmp;
}

#endif // ROTATION
#endif
