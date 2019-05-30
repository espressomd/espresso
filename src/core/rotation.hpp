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
#ifndef ROTATION_H
#define ROTATION_H
/** \file
    This file contains all subroutines required to process rotational motion.

*/

#include "config.hpp"

#ifdef ROTATION

#include "particle_data.hpp"
#include <utils/Vector.hpp>

constexpr const int ROTATION_X = 2;
constexpr const int ROTATION_Y = 4;
constexpr const int ROTATION_Z = 8;

/*************************************************************
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/

/** Propagate angular velocities and update quaternions on a particle */
void propagate_omega_quat_particle(Particle *p);

/** Convert torques to the body-fixed frame and propagate
    angular velocities */
void convert_torques_propagate_omega();

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques();

Utils::Vector3d convert_vector_body_to_space(const Particle &p,
                                             const Utils::Vector3d &v);
Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                             const Utils::Vector3d &v);

inline void convert_quat_to_director(const Utils::Vector4d &quat,
                                     Utils::Vector3d &director) {
  /* director */
  director[0] = 2 * (quat[1] * quat[3] + quat[0] * quat[2]);
  director[1] = 2 * (quat[2] * quat[3] - quat[0] * quat[1]);
  director[2] = (quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] +
                 quat[3] * quat[3]);
}

/** Multiply two quaternions */
template <typename T1, typename T2, typename T3>
void multiply_quaternions(const T1 &a, const T2 &b, T3 &result) {
  // Formula from http://www.j3d.org/matrix_faq/matrfaq_latest.html
  result[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  result[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
  result[2] = a[0] * b[2] + a[2] * b[0] + a[3] * b[1] - a[1] * b[3];
  result[3] = a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1];
}

/** Convert director to quaternions */
int convert_director_to_quat(const Utils::Vector3d &d, Utils::Vector4d &quat);

#ifdef DIPOLES

/** convert a dipole moment to quaternions and dipolar strength  */
inline std::pair<Utils::Vector4d, double>
convert_dip_to_quat(const Utils::Vector3d &dip) {
  Utils::Vector4d quat;
  convert_director_to_quat(dip, quat);

  return {quat, dip.norm()};
}

#endif

/** Rotate the particle p around the NORMALIZED axis a by amount phi */
void local_rotate_particle(Particle &p, const Utils::Vector3d &a, double phi);

inline void normalize_quaternion(double *q) {
  double tmp = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  q[0] /= tmp;
  q[1] /= tmp;
  q[2] /= tmp;
  q[3] /= tmp;
}

#ifdef BROWNIAN_DYNAMICS
/** Propagate quaternions: viscous drag driven by conservative torques.*/
void bd_drag_rot(Particle &p, double dt);
/** Set the terminal angular velocity driven by the conservative torques drag.*/
void bd_drag_vel_rot(Particle &p, double dt);

/** Propagate quaternion: random walk part.*/
void bd_random_walk_rot(Particle &p, double dt);
/** Thermalize angular velocity: random walk part.*/
void bd_random_walk_vel_rot(Particle &p, double dt);
#endif // BROWNIAN_DYNAMICS

#endif
#endif
