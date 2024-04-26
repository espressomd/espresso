/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ROTATION

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/mask.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/matrix.hpp>
#include <utils/quaternion.hpp>

#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

/** @brief Propagate angular velocities and update quaternions on a
 *  particle.
 */
void propagate_omega_quat_particle(Particle &p, double time_step);

void convert_torque_propagate_omega(Particle &p, double time_step);

/** Convert torques to the body-fixed frame before the integration loop. */
void convert_initial_torques(const ParticleRange &particles);

// Frame conversion routines
inline Utils::Vector3d
convert_vector_body_to_space(const Particle &p, const Utils::Vector3d &vec) {
  return p.quat() * vec;
}

inline Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                                    const Utils::Vector3d &v) {
  assert(p.quat().norm() > 0.0);
  return rotation_matrix(p.quat()).transposed() * v;
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
 *     A' = O A O^T.
 * \f]
 *
 * @tparam T Scalar type
 * @param quat quaternion to transform from, i.e. the rotation
 *             that transforms space- to body-fixed frame.
 * @param A Matrix representation in body-fixed coordinates.
 * @return Matrix representation in space-fixed coordinates.
 */
template <class T>
auto convert_body_to_space(const Utils::Quaternion<double> &quat,
                           const Utils::Matrix<T, 3, 3> &A) {
  auto const O = rotation_matrix(quat);
  return O * A * O.transposed();
}

/**
 * @brief Transform matrix from body- to space-fixed frame.
 * @tparam T Scalar type
 * @param p Particle transforming from.
 * @param A Matrix representation in body-fixed coordinates.
 * @return Matrix representation in space-fixed coordinates.
 */
template <class T>
auto convert_body_to_space(const Particle &p, const Utils::Matrix<T, 3, 3> &A) {
  return convert_body_to_space(p.quat(), A);
}

#ifdef DIPOLES

/** convert a dipole moment to quaternions and dipolar strength  */
inline std::pair<Utils::Quaternion<double>, double>
convert_dip_to_quat(const Utils::Vector3d &dip) {
  auto quat = Utils::convert_director_to_quaternion(dip);
  return {quat, dip.norm()};
}

#endif

/** Rotate the particle p around the body-frame defined NORMALIZED axis
 *  @p aBodyFrame by amount @p phi.
 */
inline Utils::Quaternion<double>
local_rotate_particle_body(Particle const &p,
                           const Utils::Vector3d &axis_body_frame,
                           const double phi) {
  // Rotation turned off entirely?
  if (!p.can_rotate())
    return p.quat();
  if (std::abs(phi) > std::numeric_limits<double>::epsilon())
    return p.quat() *
           boost::qvm::rot_quat(mask(p.rotation(), axis_body_frame), phi);
  return p.quat();
}

/** Rotate the particle p around the NORMALIZED axis aSpaceFrame by amount phi
 */
inline void local_rotate_particle(Particle &p,
                                  const Utils::Vector3d &axis_space_frame,
                                  const double phi) {
  if (std::abs(phi) > std::numeric_limits<double>::epsilon()) {
    // Convert rotation axis to body-fixed frame
    Utils::Vector3d axis = convert_vector_space_to_body(p, axis_space_frame);
    p.quat() = local_rotate_particle_body(p, axis, phi);
  }
}

inline void convert_torque_to_body_frame_apply_fix(Particle &p) {
  auto const torque = convert_vector_space_to_body(p, p.torque());
  p.torque() = mask(p.rotation(), torque);
}

#endif // ROTATION
#endif
