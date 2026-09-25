/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

/** \file
 *  This file contains all subroutines required to process rotational motion.
 */

#include <config/config.hpp>

#ifdef ESPRESSO_ROTATION

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/mask.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/matrix.hpp>
#include <utils/quaternion.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

/** @brief Propagate angular velocities and update quaternions (value form).
 *
 *  VV-rotation column-kernel core: operates on quaternion / omega value
 *  references (in/out) plus rinertia / torque / rotation values read from
 *  the store columns by the caller. The omega is axis-masked internally
 *  before the quaternion-derivative step; the caller writes the masked omega
 *  back to the column BEFORE this runs (see the VV-rotation kernel below).
 */
void propagate_omega_quat_values(Utils::Quaternion<double> &quat,
                                 Utils::Vector3d &omega,
                                 Utils::Vector3d const &rinertia,
                                 Utils::Vector3d const &torque,
                                 std::uint8_t rotation, double time_step);

/** @brief Convert torques and propagate angular velocities (value form).
 *
 *  VV-rotation column-kernel core. The caller has ALREADY converted the torque
 *  to the body frame and applied the fix (via
 *  @ref convert_torque_to_body_frame_apply_fix) and passes that torque here.
 */
void convert_torque_propagate_omega_values(Utils::Vector3d &omega,
                                           Utils::Vector3d const &rinertia,
                                           Utils::Vector3d const &torque,
                                           double time_step);

/** @brief Propagate angular velocities and update quaternions on a particle.
 *  View-form wrapper over @ref propagate_omega_quat_values; retained for the
 *  symplectic-Euler rotation path (still on the view path).
 */
void propagate_omega_quat_particle(Particle &p, double time_step);

/** View-form wrapper over @ref convert_torque_propagate_omega_values. */
void convert_torque_propagate_omega(Particle &p, double time_step);

/** Convert torques to the body-fixed frame before the integration loop. */
void convert_initial_torques(const ParticleRange &particles);

// Frame conversion routines
inline Utils::Vector3d
convert_vector_body_to_space(const Particle &p, const Utils::Vector3d &vec) {
  return Utils::Quaternion<double>(p.quat()) * vec;
}

inline Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                                    const Utils::Vector3d &v) {
  auto const quaternion = Utils::Quaternion<double>(p.quat());
  assert(quaternion.norm() > 0.0);
  return rotation_matrix(quaternion).transposed() * v;
}

// Value-parameter overloads of the frame-conversion routines: the Brownian
// column kernels have the quaternion in a local (read from the column once) so
// they pass it directly instead of a Particle view.
inline Utils::Vector3d
convert_vector_body_to_space(Utils::Quaternion<double> const &quat,
                             const Utils::Vector3d &vec) {
  return quat * vec;
}

inline Utils::Vector3d
convert_vector_space_to_body(Utils::Quaternion<double> const &quat,
                             const Utils::Vector3d &v) {
  assert(quat.norm() > 0.0);
  return rotation_matrix(quat).transposed() * v;
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

#ifdef ESPRESSO_DIPOLES

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
    return Utils::Quaternion<double>(p.quat());
  if (std::abs(phi) > std::numeric_limits<double>::epsilon())
    return Utils::Quaternion<double>(p.quat()) *
           boost::qvm::rot_quat(mask(p.rotation(), axis_body_frame), phi);
  return Utils::Quaternion<double>(p.quat());
}

/** @brief Value-parameter overload of @ref local_rotate_particle_body.
 *
 *  The Brownian rotation column kernel has the quaternion and the rotation
 *  bitfield in locals (read from the columns once) and passes them directly.
 *  The Particle-view overload is retained for other callers.
 */
inline Utils::Quaternion<double> local_rotate_particle_body(
    Utils::Quaternion<double> const &quat, std::uint8_t const rotation,
    const Utils::Vector3d &axis_body_frame, const double phi) {
  // Rotation turned off entirely?
  if (rotation == 0u)
    return quat;
  if (std::abs(phi) > std::numeric_limits<double>::epsilon())
    return quat * boost::qvm::rot_quat(mask(rotation, axis_body_frame), phi);
  return quat;
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
  auto torque_ref = p.torque();
  auto const torque = convert_vector_space_to_body(p, torque_ref);
  torque_ref = mask(p.rotation(), torque);
}

/** @brief Value-parameter overload of @ref
 * convert_torque_to_body_frame_apply_fix.
 *
 *  Converts the space-frame @p torque to the body frame using @p quat and
 *  applies the rotation-axis fix, returning the masked body-frame torque.
 *  The column kernel reads quat / torque / rotation from the store columns
 *  and writes the result back.
 */
inline Utils::Vector3d
convert_torque_to_body_frame_apply_fix(Utils::Quaternion<double> const &quat,
                                       Utils::Vector3d const &torque,
                                       std::uint8_t const rotation) {
  auto const torque_body = rotation_matrix(quat).transposed() * torque;
  return mask(rotation, torque_body);
}

// -- VV-rotation column kernels -------------------------------------------
// The caller (integrator_step_1 / _2) hoists the quaternion / angular-velocity
// / rinertia / torque / rotation *_view() handles ONCE and calls these per-row
// kernels. Columns are read into locals, the value-form cores run the math,
// and results are written back.

template <class QuatView, class OmegaView, class RinertiaView, class TorqueView,
          class RotationView>
inline void velocity_verlet_rotator_1_kernel(QuatView const &quat_view,
                                             OmegaView const &omega_view,
                                             RinertiaView const &rinertia_view,
                                             TorqueView const &torque_view,
                                             RotationView const &rotation_view,
                                             int const row, double time_step) {
  std::uint8_t const rotation = rotation_view(row);
  if (rotation == 0u)
    return; // matches Particle::can_rotate()

  Utils::Quaternion<double> quat;
  for (unsigned int j = 0u; j < 4u; ++j)
    quat[j] = quat_view(row, j);
  Utils::Vector3d omega{omega_view(row, 0), omega_view(row, 1),
                        omega_view(row, 2)};
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d const rinertia{rinertia_view(row, 0), rinertia_view(row, 1),
                                 rinertia_view(row, 2)};
#else
  Utils::Vector3d const rinertia{1., 1., 1.};
#endif
  Utils::Vector3d const torque{torque_view(row, 0), torque_view(row, 1),
                               torque_view(row, 2)};

  // The masked omega must be written back to the column BEFORE the quaternion
  // update. propagate_omega_quat_values masks the local omega first, so
  // writing that masked value back here establishes the correct ordering;
  // the value core then finishes the omega/quaternion propagation on the same
  // local.
  auto const omega_masked = Utils::mask(rotation, omega);
  for (unsigned int j = 0u; j < 3u; ++j)
    omega_view(row, j) = omega_masked[j];

  propagate_omega_quat_values(quat, omega, rinertia, torque, rotation,
                              time_step);

  for (unsigned int j = 0u; j < 3u; ++j)
    omega_view(row, j) = omega[j];
  for (unsigned int j = 0u; j < 4u; ++j)
    quat_view(row, j) = quat[j];
}

template <class QuatView, class OmegaView, class RinertiaView, class TorqueView,
          class RotationView>
inline void velocity_verlet_rotator_2_kernel(QuatView const &quat_view,
                                             OmegaView const &omega_view,
                                             RinertiaView const &rinertia_view,
                                             TorqueView const &torque_view,
                                             RotationView const &rotation_view,
                                             int const row, double time_step) {
  std::uint8_t const rotation = rotation_view(row);
  if (rotation == 0u)
    return; // matches Particle::can_rotate()

  Utils::Quaternion<double> quat;
  for (unsigned int j = 0u; j < 4u; ++j)
    quat[j] = quat_view(row, j);
  Utils::Vector3d const torque_space{torque_view(row, 0), torque_view(row, 1),
                                     torque_view(row, 2)};
  // Convert the torque to the body frame and apply the axis fix, writing it
  // back (matches convert_torque_to_body_frame_apply_fix(p) at the top of the
  // pre-8a convert_torque_propagate_omega).
  auto const torque =
      convert_torque_to_body_frame_apply_fix(quat, torque_space, rotation);
  for (unsigned int j = 0u; j < 3u; ++j)
    torque_view(row, j) = torque[j];

  Utils::Vector3d omega{omega_view(row, 0), omega_view(row, 1),
                        omega_view(row, 2)};
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d const rinertia{rinertia_view(row, 0), rinertia_view(row, 1),
                                 rinertia_view(row, 2)};
#else
  Utils::Vector3d const rinertia{1., 1., 1.};
#endif

  convert_torque_propagate_omega_values(omega, rinertia, torque, time_step);

  for (unsigned int j = 0u; j < 3u; ++j)
    omega_view(row, j) = omega[j];
}

#endif // ESPRESSO_ROTATION
