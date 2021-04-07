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
 *  A velocity Verlet algorithm using quaternions is implemented to tackle
 *  rotational motion. See @cite martys99a for the method and
 *  @cite allen2017 for the quaternion components indexing used here.
 *  A random torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature
 *  all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of
 *  freedom if ROTATION is compiled in.
 */

#include "rotation.hpp"

#ifdef ROTATION

#include <utils/Vector.hpp>
#include <utils/mask.hpp>

#include <cassert>
#include <cmath>

/** @brief Calculate the derivatives of the quaternion and angular
 *  acceleration for a given particle.
 *  See @cite sonnenschein85a. Please note that ESPResSo uses scalar-first
 *  notation for quaternions, while @cite sonnenschein85a uses scalar-last
 *  notation.
 *  @param[in]  p    %Particle
 *  @param[out] Qd   First derivative of the particle quaternion
 *  @param[out] Qdd  Second derivative of the particle quaternion
 *  @param[out] S    Function of @p Qd and @p Qdd, used to evaluate the
 *                   Lagrange parameter lambda
 *  @param[out] Wd   Angular acceleration of the particle
 */
static void define_Qdd(Particle const &p, Utils::Quaternion<double> &Qd,
                       Utils::Quaternion<double> &Qdd, Utils::Vector3d &S,
                       Utils::Vector3d &Wd) {
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
  if (p.p.rotation & ROTATION_Y)
    Wd[1] = (p.f.torque[1] + p.m.omega[2] * p.m.omega[0] *
                                 (p.p.rinertia[2] - p.p.rinertia[0])) /
            p.p.rinertia[1];
  if (p.p.rotation & ROTATION_Z)
    Wd[2] = (p.f.torque[2] + p.m.omega[0] * p.m.omega[1] *
                                 (p.p.rinertia[0] - p.p.rinertia[1])) /
            p.p.rinertia[2];

  auto const S1 = Qd.norm2();

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
  S[1] = Utils::dot(Qd, Qdd);
  S[2] = Qdd.norm2();
}

/**
 *  See @cite omelyan98a. Please note that ESPResSo uses scalar-first
 *  notation for quaternions, while @cite omelyan98a uses scalar-last
 *  notation.
 *
 *  For very high angular velocities (e.g. if the product of @p time_step
 *  with the largest component of @ref ParticleMomentum::omega "p.m.omega"
 *  is superior to ~2.0), the calculation might fail.
 *
 *  \todo implement for fixed_coord_flag
 */
void propagate_omega_quat_particle(Particle &p, double time_step) {

  // If rotation for the particle is disabled entirely, return early.
  if (p.p.rotation == ROTATION_FIXED)
    return;

  Utils::Quaternion<double> Qd{}, Qdd{};
  Utils::Vector3d S{}, Wd{};

  // Clear rotational velocity for blocked rotation axes.
  p.m.omega = Utils::mask(p.p.rotation, p.m.omega);

  define_Qdd(p, Qd, Qdd, S, Wd);

  auto const time_step_squared = time_step * time_step;
  auto const time_step_half = 0.5 * time_step;

  /* Eq. (12) @cite omelyan98a. */
  auto const square =
      1 - time_step_squared *
              (S[0] +
               time_step * (S[1] + time_step_half / 2. * (S[2] - S[0] * S[0])));
  assert(square >= 0.);
  auto const lambda = 1 - S[0] * 0.5 * time_step_squared - sqrt(square);

  p.m.omega += time_step_half * Wd;
  p.r.quat += time_step * (Qd + time_step_half * Qdd) - lambda * p.r.quat;

  /* and rescale quaternion, so it is exactly of unit length */
  auto const scale = p.r.quat.norm();
  if (scale == 0) {
    p.r.quat = Utils::Quaternion<double>::identity();
  } else {
    p.r.quat /= scale;
  }
}

void convert_torques_propagate_omega(const ParticleRange &particles,
                                     double time_step) {
  for (auto &p : particles) {
    // Skip particle if rotation is turned off entirely for it.
    if (p.p.rotation == ROTATION_FIXED)
      continue;

    convert_torque_to_body_frame_apply_fix(p);

    // Propagation of angular velocities
    p.m.omega += hadamard_division(0.5 * time_step * p.f.torque, p.p.rinertia);

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

      p.m.omega = omega_0 + (0.5 * time_step) * Wd;
    }
  }
}

void convert_initial_torques(const ParticleRange &particles) {
  for (auto &p : particles) {
    if (!p.p.rotation)
      continue;
    convert_torque_to_body_frame_apply_fix(p);
  }
}

#endif // ROTATION
