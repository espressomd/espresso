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
 *  rotational motion. See "Velocity Verlet algorithm for
 *  dissipative-particle-dynamics-based models of suspensions", Martys
 *  and Mountain 1999 (10.1103/PhysRevE.59.3733) for the method and
 *  "Computer simulation of liquids", Allen and Tildesley 2017
 *  (10.1093/oso/9780198803195.001.0001) for the quaternion components
 *  indexing used here.
 *  A random torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature
 * all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of
 * freedom if ROTATION
 *  flag is set in \ref config.hpp "config.hpp".
 */

#include "rotation.hpp"

#ifdef ROTATION
#include "cells.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "forces.hpp"
#include "global.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"

#include <utils/constants.hpp>
#include <utils/math/rotation_matrix.hpp>

#include <cmath>

/** Calculate the derivatives of the quaternion and angular acceleration
 *  for a given particle.
 *  See "An improved algorithm for molecular dynamics simulation of rigid
 *  molecules", Sonnenschein and Roland 1985 (10.1016/0021-9991(85)90151-2).
 *  Please note that ESPResSo uses scalar-first notation for quaternions,
 *  while the paper uses scalar-last notation.
 *  @param[in]  p    %Particle
 *  @param[out] Qd   First derivative of the particle quaternion
 *  @param[out] Qdd  Second derivative of the particle quaternion
 *  @param[out] S    Function of @p Qd and @p Qdd, used to evaluate the
 *                   Lagrange parameter lambda
 *  @param[out] Wd   Angular acceleration of the particle
 */
static void define_Qdd(Particle const &p, Utils::Vector4d &Qd,
                       Utils::Vector4d &Qdd, Utils::Vector3d &S,
                       Utils::Vector3d &Wd);

/** Convert director to quaternions */
int convert_director_to_quat(const Utils::Vector3d &d, Utils::Vector4d &quat) {
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

void define_Qdd(Particle const &p, Utils::Vector4d &Qd, Utils::Vector4d &Qdd,
                Utils::Vector3d &S, Utils::Vector3d &Wd) {
  /* calculate the first derivative of the quaternion */
  /* Taken from "An improved algorithm for molecular dynamics simulation of
   * rigid molecules", Sonnenschein, Roland (1985), Eq. 4.*/
  Qd[0] = 0.5 * (-p.r.quat[1] * p.m.omega[0] - p.r.quat[2] * p.m.omega[1] -
                 p.r.quat[3] * p.m.omega[2]);

  Qd[1] = 0.5 * (p.r.quat[0] * p.m.omega[0] - p.r.quat[3] * p.m.omega[1] +
                 p.r.quat[2] * p.m.omega[2]);

  Qd[2] = 0.5 * (p.r.quat[3] * p.m.omega[0] + p.r.quat[0] * p.m.omega[1] -
                 p.r.quat[1] * p.m.omega[2]);

  Qd[3] = 0.5 * (-p.r.quat[2] * p.m.omega[0] + p.r.quat[1] * p.m.omega[1] +
                 p.r.quat[0] * p.m.omega[2]);

  /* Calculate the angular acceleration. */
  /* Taken from "An improved algorithm for molecular dynamics simulation of
   * rigid molecules", Sonnenschein, Roland (1985), Eq. 5.*/
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
  /* Taken from "An improved algorithm for molecular dynamics simulation of
   * rigid molecules", Sonnenschein, Roland (1985), Eq. 8.*/
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
  S[1] = Qd * Qdd;
  S[2] = Qdd.norm2();
}

/** Propagate angular velocities and quaternions.
 *  See "On the numerical integration of motion for rigid polyatomics:
 *  The modified quaternion approach", Omelyan 1998 (10.1063/1.168642).
 *  Please note that ESPResSo uses scalar-first notation for quaternions,
 *  while the paper uses scalar-last notation.
 *
 *  For very high angular velocities (e.g. if the product of @ref time_step
 *  with the largest component of @ref ParticleMomentum::omega "p.m.omega"
 *  is superior to ~2.0), the calculation might fail.
 * \todo implement for fixed_coord_flag
 */
void propagate_omega_quat_particle(Particle &p) {

  // If rotation for the particle is disabled entirely, return early.
  if (!p.p.rotation)
    return;

  Utils::Vector4d Qd{}, Qdd{};
  Utils::Vector3d S{}, Wd{};

  // Clear rotational velocity for blocked rotation axes.
  if (!(p.p.rotation & ROTATION_X))
    p.m.omega[0] = 0;
  if (!(p.p.rotation & ROTATION_Y))
    p.m.omega[1] = 0;
  if (!(p.p.rotation & ROTATION_Z))
    p.m.omega[2] = 0;

  define_Qdd(p, Qd, Qdd, S, Wd);

  /* Taken from "On the numerical integration of motion for rigid polyatomics:
   * The modified quaternion approach", Omelyan (1998), Eq. 12.*/
  auto const square =
      1 - time_step_squared *
              (S[0] +
               time_step * (S[1] + time_step_half / 2. * (S[2] - S[0] * S[0])));
  assert(square >= 0.);
  auto const lambda = 1 - S[0] * time_step_squared_half - sqrt(square);

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

inline void convert_torque_to_body_frame_apply_fix_and_thermostat(Particle &p) {
  auto const t = convert_vector_space_to_body(p, p.f.torque);
  p.f.torque = Utils::Vector3d{};

  if (thermo_switch & THERMO_LANGEVIN) {
#if defined(VIRTUAL_SITES) && defined(THERMOSTAT_IGNORE_NON_VIRTUAL)
    if (!p.p.is_virtual)
#endif
    {
      friction_thermo_langevin_rotation(p);

      p.f.torque += t;
    }
  } else {
    p.f.torque = t;
  }

  if (!(p.p.rotation & ROTATION_X))
    p.f.torque[0] = 0;

  if (!(p.p.rotation & ROTATION_Y))
    p.f.torque[1] = 0;

  if (!(p.p.rotation & ROTATION_Z))
    p.f.torque[2] = 0;
}

/** convert the torques to the body-fixed frames and propagate angular
 * velocities */
void convert_torques_propagate_omega(const ParticleRange &particles) {

#if defined(CUDA) && defined(ENGINE)
  if ((lb_lbfluid_get_lattice_switch() == ActiveLB::GPU) &&
      swimming_particles_exist) {
    copy_v_cs_from_GPU(particles);
  }
#endif

  for (auto &p : particles) {
    // Skip particle if rotation is turned off entirely for it.
    if (!p.p.rotation)
      continue;

    convert_torque_to_body_frame_apply_fix_and_thermostat(p);

#if defined(ENGINE)
    if (p.swim.swimming && lb_lbfluid_get_lattice_switch() != ActiveLB::NONE) {
      if (lb_lbfluid_get_lattice_switch() == ActiveLB::CPU && n_nodes > 1 &&
          p.swim.rotational_friction != 0.) {
        runtimeErrorMsg() << "ENGINE rotational_friction feature with CPU-LB "
                             "only implemented for one CPU core";
      }

      auto const dip = p.swim.dipole_length * p.r.calc_director();

      auto const diff = p.swim.v_center - p.swim.v_source;

      const Utils::Vector3d cross = vector_product(diff, dip);
      const double l_diff = diff.norm();
      const double l_cross = cross.norm();

      if (l_cross > 0 && p.swim.dipole_length > 0) {
        auto const omega_swim =
            l_diff / (l_cross * p.swim.dipole_length) * cross;

        auto const omega_swim_body =
            convert_vector_space_to_body(p, omega_swim);
        p.f.torque +=
            p.swim.rotational_friction * (omega_swim_body - p.m.omega);
      }
    }
#endif

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
    convert_torque_to_body_frame_apply_fix_and_thermostat(p);
  }
}
// Frame conversion routines

Utils::Vector3d convert_vector_body_to_space(const Particle &p,
                                             const Utils::Vector3d &vec) {
  auto const A = rotation_matrix(p.r.quat);
  return transpose(A) * vec;
}

Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                             const Utils::Vector3d &v) {
  return rotation_matrix(p.r.quat) * v;
}

/** Rotate the particle p around the NORMALIZED axis aSpaceFrame by amount phi
 */
void local_rotate_particle(Particle &p, const Utils::Vector3d &axis_space_frame,
                           const double phi) {
  // Convert rotation axis to body-fixed frame
  Utils::Vector3d axis = convert_vector_space_to_body(p, axis_space_frame);

  // Rotation turned off entirely?
  if (!p.p.rotation)
    return;

  // Per coordinate fixing
  if (!(p.p.rotation & ROTATION_X))
    axis[0] = 0;
  if (!(p.p.rotation & ROTATION_Y))
    axis[1] = 0;
  if (!(p.p.rotation & ROTATION_Z))
    axis[2] = 0;
  // Re-normalize rotation axis
  double l = axis.norm();
  // Check, if the rotation axis is nonzero
  if (l < std::numeric_limits<double>::epsilon())
    return;

  axis /= l;

  double q[4];
  q[0] = cos(phi / 2);
  double tmp = sin(phi / 2);
  q[1] = tmp * axis[0];
  q[2] = tmp * axis[1];
  q[3] = tmp * axis[2];

  // Normalize
  normalize_quaternion(q);

  // Rotate the particle
  double qn[4]; // Resulting quaternion
  multiply_quaternions(p.r.quat, q, qn);
  for (int k = 0; k < 4; k++)
    p.r.quat[k] = qn[k];
}

#endif
