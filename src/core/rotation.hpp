/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef ROTATION_H
#define ROTATION_H
/** \file
 *  This file contains all subroutines required to process rotational motion.
 */

#include "config.hpp"

#ifdef ROTATION

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/mask.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/rotation_matrix.hpp>

/** Propagate angular velocities and update quaternions on a particle */
void propagate_omega_quat_particle(Particle &p);

/** Convert torques to the body-fixed frame and propagate
    angular velocities */
void convert_torques_propagate_omega(const ParticleRange &particles);

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques(const ParticleRange &particles);

// Frame conversion routines
inline Utils::Vector3d
convert_vector_body_to_space(const Particle &p, const Utils::Vector3d &vec) {
  return rotation_matrix(p.r.quat) * vec;
}

inline Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                                    const Utils::Vector3d &v) {
  return transpose(rotation_matrix(p.r.quat)) * v;
}

/**
 * @brief Transform matrix from body- to space-fixed frame.
 *
 * Given a linear map represented by \f$ A \in \mathbb{R}^{3 \times 3}\f$
 * in the body-fixed frame, this returns the matrix \f$ A \in \mathbb{R}^{3
 * \times 3}\f$ representing the map in the space-fixed frame. They are related
 * by the map between the space-fixed and body-fixed frame \f$O\f$ like
 *
 * \f[
 *     A' = O^T A O.
 * \f]
 *
 * @tparam T Scalar type
 * @param p Particle transforming from.
 * @param A Matrix to transform
 * @return Matrix representation in space-fixed coordinates.
 */
template <class T>
auto convert_body_to_space(const Particle &p, const Utils::Matrix<T, 3, 3> &A) {
  auto const O = rotation_matrix(p.r.quat);
  return transpose(O) * A * O;
}

#ifdef DIPOLES

/** convert a dipole moment to quaternions and dipolar strength  */
inline std::pair<Utils::Vector4d, double>
convert_dip_to_quat(const Utils::Vector3d &dip) {
  auto quat = Utils::convert_director_to_quaternion(dip);
  return {quat, dip.norm()};
}

#endif

/** Rotate the particle p around the body-frame defined NORMALIZED axis
 *  @p aBodyFrame by amount @p phi.
 */
inline Utils::Vector4d
local_rotate_particle_body(Particle const &p,
                           const Utils::Vector3d &axis_body_frame,
                           const double phi) {
  auto axis = axis_body_frame;

  // Rotation turned off entirely?
  if (!p.p.rotation)
    return {};

  // Convert rotation axis to body-fixed frame
  axis = mask(p.p.rotation, axis).normalize();

  auto const s = std::sin(phi / 2);
  auto const q =
      Utils::Vector4d{cos(phi / 2), s * axis[0], s * axis[1], s * axis[2]}
          .normalize();

  return Utils::multiply_quaternions(p.r.quat, q);
}

/** Rotate the particle p around the NORMALIZED axis aSpaceFrame by amount phi
 */
inline void local_rotate_particle(Particle &p,
                                  const Utils::Vector3d &axis_space_frame,
                                  const double phi) {
  // Convert rotation axis to body-fixed frame
  Utils::Vector3d axis = convert_vector_space_to_body(p, axis_space_frame);
  p.r.quat = local_rotate_particle_body(p, axis, phi);
}

inline void convert_torque_to_body_frame_apply_fix(Particle &p) {
  auto const torque = convert_vector_space_to_body(p, p.f.torque);
  p.f.torque = mask(p.p.rotation, torque);
}

#endif // ROTATION
#endif
