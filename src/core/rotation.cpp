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
/** \file
 *  Molecular dynamics integrator for rotational motion.
 *
 *  A velocity Verlet <a
 * HREF="http://ciks.cbt.nist.gov/~garbocz/dpd1/dpd.html">algorithm</a>
 *  using quaternions is implemented to tackle rotational motion.
 *  See @cite allen2017 for the quaternion components indexing used here.
 *  A random torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature
 * all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of
 * freedom if ROTATION
 *  flag is set in \ref config.hpp "config.hpp".
 */

#include "rotation.hpp"

#ifdef ROTATION
#include "integrate.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/mask.hpp>

#include <cmath>

/** Calculate the derivatives of the quaternion and angular acceleration
 *  for a given particle.
 *  See @cite sonnenschein85a
 *  @param[in]  p    %Particle
 *  @param[out] Qd   First derivative of the particle quaternion
 *  @param[out] Qdd  Second derivative of the particle quaternion
 *  @param[out] S    Function of @p Qd and @p Qdd, used to evaluate the
 *                   Lagrange parameter lambda
 *  @param[out] Wd   Angular acceleration of the particle
 */
static void define_Qdd(Particle const &p, double Qd[4], double Qdd[4],
                       double S[3], double Wd[3]) {
  /* calculate the first derivative of the quaternion */
  /* Eq. (4) @cite sonnenschein85a */
  Qd[0] = 0.5 * (-p.r.quat[1] * p.m.omega[0] - p.r.quat[2] * p.m.omega[1] -
                 p.r.quat[3] * p.m.omega[2]);

  Qd[1] = 0.5 * (p.r.quat[0] * p.m.omega[0] - p.r.quat[3] * p.m.omega[1] +
                 p.r.quat[2] * p.m.omega[2]);

  Qd[2] = 0.5 * (p.r.quat[3] * p.m.omega[0] + p.r.quat[0] * p.m.omega[1] -
                 p.r.quat[1] * p.m.omega[2]);

  Qd[3] = 0.5 * (-p.r.quat[2] * p.m.omega[0] + p.r.quat[1] * p.m.omega[1] +
                 p.r.quat[0] * p.m.omega[2]);

  /* Calculate the angular acceleration. */
  /* Eq. (5) @cite sonnenschein85a */
  if (p.p.rotation & ROTATION_X)
    Wd[0] = (p.f.torque[0] + p.m.omega[1] * p.m.omega[2] *
                                 (p.p.rinertia[1] - p.p.rinertia[2])) /
            p.p.rinertia[0];
  else
    Wd[0] = 0.0;
  if (p.p.rotation & ROTATION_Y)
    Wd[1] = (p.f.torque[1] + p.m.omega[2] * p.m.omega[0] *
                                 (p.p.rinertia[2] - p.p.rinertia[0])) /
            p.p.rinertia[1];
  else
    Wd[1] = 0.0;
  if (p.p.rotation & ROTATION_Z)
    Wd[2] = (p.f.torque[2] + p.m.omega[0] * p.m.omega[1] *
                                 (p.p.rinertia[0] - p.p.rinertia[1])) /
            p.p.rinertia[2];
  else
    Wd[2] = 0.0;

  auto const S1 = Qd[0] * Qd[0] + Qd[1] * Qd[1] + Qd[2] * Qd[2] + Qd[3] * Qd[3];

  /* Calculate the second derivative of the quaternion. */
  /* Eq. (8) @cite sonnenschein85a */
  Qdd[0] =
      0.5 * (-p.r.quat[1] * Wd[0] - p.r.quat[2] * Wd[1] - p.r.quat[3] * Wd[2]) -
      p.r.quat[0] * S1;

  Qdd[1] =
      0.5 * (p.r.quat[0] * Wd[0] - p.r.quat[3] * Wd[1] + p.r.quat[2] * Wd[2]) -
      p.r.quat[1] * S1;

  Qdd[2] =
      0.5 * (p.r.quat[3] * Wd[0] + p.r.quat[0] * Wd[1] - p.r.quat[1] * Wd[2]) -
      p.r.quat[2] * S1;

  Qdd[3] =
      0.5 * (-p.r.quat[2] * Wd[0] + p.r.quat[1] * Wd[1] + p.r.quat[0] * Wd[2]) -
      p.r.quat[3] * S1;

  S[0] = S1;
  S[1] = Qd[0] * Qdd[0] + Qd[1] * Qdd[1] + Qd[2] * Qdd[2] + Qd[3] * Qdd[3];
  S[2] = Qdd[0] * Qdd[0] + Qdd[1] * Qdd[1] + Qdd[2] * Qdd[2] + Qdd[3] * Qdd[3];
}

/** propagate angular velocities and quaternions
 * \todo implement for fixed_coord_flag
 */
void propagate_omega_quat_particle(Particle &p) {

  // If rotation for the particle is disabled entirely, return early.
  if (p.p.rotation == ROTATION_FIXED)
    return;

  Utils::Vector4d Qd{}, Qdd{};
  Utils::Vector3d S{}, Wd{};

  // Clear rotational velocity for blocked rotation axes.
  p.m.omega = Utils::mask(p.p.rotation, p.m.omega);

  define_Qdd(p, Qd.data(), Qdd.data(), S.data(), Wd.data());

  /* Eq. (12) @cite sonnenschein85a. */
  auto const lambda =
      1 - S[0] * time_step_squared_half -
      sqrt(1 - time_step_squared *
                   (S[0] + time_step * (S[1] + time_step_half / 2. *
                                                   (S[2] - S[0] * S[0]))));

  p.m.omega += time_step_half * Wd;
  p.r.quat += time_step * (Qd + time_step_half * Qdd) - lambda * p.r.quat;

  /* and rescale quaternion, so it is exactly of unit length */
  auto const scale = p.r.quat.norm();
  if (scale == 0) {
    p.r.quat[0] = 1;
  } else {
    p.r.quat /= scale;
  }
}

/** convert the torques to the body-fixed frames and propagate angular
 * velocities */
void convert_torques_propagate_omega(const ParticleRange &particles) {
  for (auto &p : particles) {
    // Skip particle if rotation is turned off entirely for it.
    if (!p.p.rotation)
      continue;

    convert_torque_to_body_frame_apply_fix(p);

    // Propagation of angular velocities
    p.m.omega += hadamard_division(time_step_half * p.f.torque, p.p.rinertia);

    // zeroth estimate of omega
    Utils::Vector3d omega_0 = p.m.omega;

    /* if the tensor of inertia is isotropic, the following refinement is not
       needed.
       Otherwise repeat this loop 2-3 times depending on the required accuracy
       */

    const double rinertia_diff_01 = p.p.rinertia[0] - p.p.rinertia[1];
    const double rinertia_diff_12 = p.p.rinertia[1] - p.p.rinertia[2];
    const double rinertia_diff_20 = p.p.rinertia[2] - p.p.rinertia[0];
    for (int times = 0; times <= 5; times++) {
      Utils::Vector3d Wd;

      Wd[0] = p.m.omega[1] * p.m.omega[2] * rinertia_diff_12 / p.p.rinertia[0];
      Wd[1] = p.m.omega[2] * p.m.omega[0] * rinertia_diff_20 / p.p.rinertia[1];
      Wd[2] = p.m.omega[0] * p.m.omega[1] * rinertia_diff_01 / p.p.rinertia[2];

      p.m.omega = omega_0 + time_step_half * Wd;
    }
  }
}

/** convert the torques to the body-fixed frames before the integration loop */
void convert_initial_torques(const ParticleRange &particles) {
  for (auto &p : particles) {
    if (!p.p.rotation)
      continue;
    convert_torque_to_body_frame_apply_fix(p);
  }
}

#endif // ROTATION
