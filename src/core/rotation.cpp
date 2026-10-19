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
/** \file
 *  Molecular dynamics integrator for rotational motion.
 *
 *  A velocity Verlet algorithm using quaternions is implemented to tackle
 *  rotational motion. See @cite martys99a for the method and
 *  @cite allen17a for the quaternion components indexing used here.
 *  A random torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature
 *  all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of
 *  freedom if ROTATION is compiled in.
 */

#include <config/config.hpp>

#ifdef ESPRESSO_ROTATION

#include "rotation.hpp"

#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/mask.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>

/** @brief Calculate the derivatives of the quaternion and angular
 *  acceleration for a given particle.
 *  See @cite sonnenschein85a. Please note that ESPResSo uses scalar-first
 *  notation for quaternions, while @cite sonnenschein85a uses scalar-last
 *  notation.
 *
 *  Value-parameter form. The VV-rotation column kernel reads the
 *  quaternion / omega / rinertia / torque / rotation columns into locals once
 *  and passes them here. The arithmetic and the per-axis @c can_rotate_around
 *  branch (a bit-test on the rotation byte) operate on the passed values.
 *  @param[in]  quaternion   Particle quaternion
 *  @param[in]  omega        Particle angular velocity (already axis-masked)
 *  @param[in]  rinertia     Particle rotational inertia
 *  @param[in]  torque       Particle torque
 *  @param[in]  rotation     Particle rotation bitfield
 *  @param[out] Qd   First derivative of the particle quaternion
 *  @param[out] Qdd  Second derivative of the particle quaternion
 *  @param[out] S    Function of @p Qd and @p Qdd, used to evaluate the
 *                   Lagrange parameter lambda
 *  @param[out] Wd   Angular acceleration of the particle
 */
static void
define_Qdd(Utils::Quaternion<double> const &quaternion,
           Utils::Vector3d const &omega, Utils::Vector3d const &rinertia,
           Utils::Vector3d const &torque, std::uint8_t const rotation,
           Utils::Quaternion<double> &Qd, Utils::Quaternion<double> &Qdd,
           Utils::Vector3d &S, Utils::Vector3d &Wd) {
  /* calculate the first derivative of the quaternion */
  /* Eq. (4) @cite sonnenschein85a */
  Qd[0] = 0.5 * (-quaternion[1] * omega[0] - quaternion[2] * omega[1] -
                 quaternion[3] * omega[2]);

  Qd[1] = 0.5 * (quaternion[0] * omega[0] - quaternion[3] * omega[1] +
                 quaternion[2] * omega[2]);

  Qd[2] = 0.5 * (quaternion[3] * omega[0] + quaternion[0] * omega[1] -
                 quaternion[1] * omega[2]);

  Qd[3] = 0.5 * (-quaternion[2] * omega[0] + quaternion[1] * omega[1] +
                 quaternion[0] * omega[2]);

  /* Calculate the angular acceleration. */
  /* Eq. (5) @cite sonnenschein85a */
  if (detail::get_nth_bit(rotation, 0))
    Wd[0] = (torque[0] + omega[1] * omega[2] * (rinertia[1] - rinertia[2])) /
            rinertia[0];
  if (detail::get_nth_bit(rotation, 1))
    Wd[1] = (torque[1] + omega[2] * omega[0] * (rinertia[2] - rinertia[0])) /
            rinertia[1];
  if (detail::get_nth_bit(rotation, 2))
    Wd[2] = (torque[2] + omega[0] * omega[1] * (rinertia[0] - rinertia[1])) /
            rinertia[2];

  auto const S1 = Qd.norm2();

  /* Calculate the second derivative of the quaternion. */
  /* Eq. (8) @cite sonnenschein85a */
  Qdd[0] = 0.5 * (-quaternion[1] * Wd[0] - quaternion[2] * Wd[1] -
                  quaternion[3] * Wd[2]) -
           quaternion[0] * S1;

  Qdd[1] = 0.5 * (quaternion[0] * Wd[0] - quaternion[3] * Wd[1] +
                  quaternion[2] * Wd[2]) -
           quaternion[1] * S1;

  Qdd[2] = 0.5 * (quaternion[3] * Wd[0] + quaternion[0] * Wd[1] -
                  quaternion[1] * Wd[2]) -
           quaternion[2] * S1;

  Qdd[3] = 0.5 * (-quaternion[2] * Wd[0] + quaternion[1] * Wd[1] +
                  quaternion[0] * Wd[2]) -
           quaternion[3] * S1;

  S[0] = S1;
  S[1] = Utils::dot(Qd, Qdd);
  S[2] = Qdd.norm2();
}

/**
 *  See @cite omelyan98a. Please note that ESPResSo uses scalar-first
 *  notation for quaternions, while @cite omelyan98a uses scalar-last
 *  notation.
 *
 *  For very high angular velocities (e.g. if the product of @p time_step
 *  with the largest component of @c omega (accessible via @ref
 * Particle::omega()) is superior to ~2.0) and for @p time_step superior or
 * equal to unity, the calculation might fail.
 *
 *  \todo implement for fixed_coord_flag
 */
void propagate_omega_quat_values(Utils::Quaternion<double> &quat,
                                 Utils::Vector3d &omega,
                                 Utils::Vector3d const &rinertia,
                                 Utils::Vector3d const &torque,
                                 std::uint8_t const rotation,
                                 double time_step) {
  assert(rotation != 0u);

  Utils::Quaternion<double> Qd{}, Qdd{};
  Utils::Vector3d S{}, Wd{};

  // Clear rotational velocity for blocked rotation axes. The masked omega MUST
  // be the value define_Qdd sees (particles with a blocked axis carrying
  // non-zero omega must see the zeroed component during quaternion-derivative
  // computation) -- the caller writes this masked value back to the omega
  // column before define_Qdd runs, reproduced here by masking the local omega
  // before the define_Qdd call below.
  omega = Utils::mask(rotation, omega);

  define_Qdd(quat, omega, rinertia, torque, rotation, Qd, Qdd, S, Wd);

  auto const time_step_squared = time_step * time_step;
  auto const time_step_half = 0.5 * time_step;

  /* Eq. (12) @cite omelyan98a. */
  auto const square =
      1 - time_step_squared *
              (S[0] +
               time_step * (S[1] + time_step_half / 2. * (S[2] - S[0] * S[0])));
  assert(square >= 0.);
  auto const lambda = 1 - S[0] * 0.5 * time_step_squared - sqrt(square);

  omega += time_step_half * Wd;
  auto const quaternion_old = quat;
  quat += time_step * (Qd + time_step_half * Qdd) - lambda * quaternion_old;

  /* and rescale quaternion, so it is exactly of unit length */
  auto const scale = quat.norm();
  if (scale == 0) {
    quat = Utils::Quaternion<double>::identity();
  } else {
    quat /= scale;
  }
}

void convert_torque_propagate_omega_values(Utils::Vector3d &omega,
                                           Utils::Vector3d const &rinertia,
                                           Utils::Vector3d const &torque,
                                           double time_step) {
  omega += hadamard_division(0.5 * time_step * torque, rinertia);

  // zeroth estimate of omega
  Utils::Vector3d omega_0 = omega;

  /* if the tensor of inertia is isotropic, the following refinement is not
     needed.
     Otherwise repeat this loop 2-3 times depending on the required accuracy
   */

  auto const rinertia_diff_01 = rinertia[0] - rinertia[1];
  auto const rinertia_diff_12 = rinertia[1] - rinertia[2];
  auto const rinertia_diff_20 = rinertia[2] - rinertia[0];
  for (int times = 0; times <= 5; times++) {
    Utils::Vector3d Wd;

    Wd[0] = omega[1] * omega[2] * rinertia_diff_12 / rinertia[0];
    Wd[1] = omega[2] * omega[0] * rinertia_diff_20 / rinertia[1];
    Wd[2] = omega[0] * omega[1] * rinertia_diff_01 / rinertia[2];

    omega = omega_0 + (0.5 * time_step) * Wd;
  }
}

// View-form wrappers over the value cores above. Retained for the non-column
// callers (symplectic-Euler rotation, which stays on the view path). They read
// the columns via the Particle view, run the shared value core, and write the
// results back -- bitwise identical to the pre-8a bodies.
void propagate_omega_quat_particle(Particle &p, double time_step) {
  assert(p.can_rotate());
  Utils::Quaternion<double> quat = p.quat();
  Utils::Vector3d omega = p.omega();
  Utils::Vector3d const rinertia = p.rinertia();
  Utils::Vector3d const torque = p.torque();
  auto const rotation = p.rotation();
  // Masked omega write-back BEFORE the quaternion-derivative step (the omega
  // column must carry the masked value when define_Qdd reads it in the pre-8a
  // body); the value core re-masks the local, so the arithmetic matches.
  p.omega() = Utils::mask(rotation, omega);
  propagate_omega_quat_values(quat, omega, rinertia, torque, rotation,
                              time_step);
  p.omega() = omega;
  p.quat() = quat;
}

void convert_torque_propagate_omega(Particle &p, double time_step) {
  assert(p.can_rotate());
  convert_torque_to_body_frame_apply_fix(p);
  Utils::Vector3d omega = p.omega();
  Utils::Vector3d const rinertia = p.rinertia();
  Utils::Vector3d const torque = p.torque();
  convert_torque_propagate_omega_values(omega, rinertia, torque, time_step);
  p.omega() = omega;
}

void convert_initial_torques(const ParticleRange &particles) {
  for (auto &p : particles) {
    if (!p.can_rotate())
      continue;
    convert_torque_to_body_frame_apply_fix(p);
  }
}

#endif // ESPRESSO_ROTATION
