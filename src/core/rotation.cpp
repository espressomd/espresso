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
/** \file
 *  Molecular dynamics integrator for rotational motion.
 *
 *  A velocity Verlet <a
 * HREF="http://ciks.cbt.nist.gov/~garbocz/dpd1/dpd.html">algorithm</a>
 *  using quaternions is implemented to tackle rotational motion. A random
 * torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature
 * all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of
 * freedom if ROTATION
 *  flag is set in \ref config.hpp "config.hpp".
 */

#include "rotation.hpp"

/****************************************************
 *                     DEFINES
 ***************************************************/
/**************** local variables  *******************/

#ifdef ROTATION
#include "brownian_inline.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "ghosts.hpp"
#include "global.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"

#include <utils/constants.hpp>
#include <utils/Vector.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

/** \name Private Functions */
/************************************************************/
/*@{*/

/** define first and second time derivatives of a quaternion, as well
    as the angular acceleration. */
static void define_Qdd(Particle *p, double Qd[4], double Qdd[4], double S[3],
                       double Wd[3]);

/*@}*/

/** convert quaternions to the director */
/** Convert director to quaternions */
int convert_director_to_quat(const Utils::Vector3d &d, Utils::Vector4d &quat) {
  double d_xy, dm;
  double theta2, phi2;

  // Calculate magnitude of the given vector
  dm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);

  // The vector needs to be != 0 to be converted into a quaternion
  if (dm < ROUND_ERROR_PREC) {
    return 1;
  }
  // Calculate angles
  d_xy = sqrt(d[0] * d[0] + d[1] * d[1]);
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
  quat[0] = cos(theta2) * cos(phi2);
  quat[1] = -sin(theta2) * cos(phi2);
  quat[2] = -sin(theta2) * sin(phi2);
  quat[3] = cos(theta2) * sin(phi2);

  return 0;
}

/** Here we use quaternions to calculate the rotation matrix which
    will be used then to transform torques from the laboratory to
    the body-fixed frames.
    Taken from "Goldstein - Classical Mechanics" (Chapter 4.6 Eq. 4.47).
*/
void define_rotation_matrix(Particle const &p, double A[9]) {
  double q0q0 = p.r.quat[0];
  q0q0 *= q0q0;

  double q1q1 = p.r.quat[1];
  q1q1 *= q1q1;

  double q2q2 = p.r.quat[2];
  q2q2 *= q2q2;

  double q3q3 = p.r.quat[3];
  q3q3 *= q3q3;

  A[0 + 3 * 0] = q0q0 + q1q1 - q2q2 - q3q3;
  A[1 + 3 * 1] = q0q0 - q1q1 + q2q2 - q3q3;
  A[2 + 3 * 2] = q0q0 - q1q1 - q2q2 + q3q3;

  A[0 + 3 * 1] = 2 * (p.r.quat[1] * p.r.quat[2] + p.r.quat[0] * p.r.quat[3]);
  A[0 + 3 * 2] = 2 * (p.r.quat[1] * p.r.quat[3] - p.r.quat[0] * p.r.quat[2]);
  A[1 + 3 * 0] = 2 * (p.r.quat[1] * p.r.quat[2] - p.r.quat[0] * p.r.quat[3]);

  A[1 + 3 * 2] = 2 * (p.r.quat[2] * p.r.quat[3] + p.r.quat[0] * p.r.quat[1]);
  A[2 + 3 * 0] = 2 * (p.r.quat[1] * p.r.quat[3] + p.r.quat[0] * p.r.quat[2]);
  A[2 + 3 * 1] = 2 * (p.r.quat[2] * p.r.quat[3] - p.r.quat[0] * p.r.quat[1]);
}

/** calculate the second derivative of the quaternion of a given particle
    as well as Wd vector which is the angular acceleration of this particle */
void define_Qdd(Particle *p, double Qd[4], double Qdd[4], double S[3],
                double Wd[3]) {
  /* calculate the first derivative of the quaternion */
  /* Taken from "An improved algorithm for molecular dynamics simulation of
   * rigid molecules", Sonnenschein, Roland (1985), Eq. 4.*/
  Qd[0] = 0.5 * (-p->r.quat[1] * p->m.omega[0] - p->r.quat[2] * p->m.omega[1] -
                 p->r.quat[3] * p->m.omega[2]);

  Qd[1] = 0.5 * (p->r.quat[0] * p->m.omega[0] - p->r.quat[3] * p->m.omega[1] +
                 p->r.quat[2] * p->m.omega[2]);

  Qd[2] = 0.5 * (p->r.quat[3] * p->m.omega[0] + p->r.quat[0] * p->m.omega[1] -
                 p->r.quat[1] * p->m.omega[2]);

  Qd[3] = 0.5 * (-p->r.quat[2] * p->m.omega[0] + p->r.quat[1] * p->m.omega[1] +
                 p->r.quat[0] * p->m.omega[2]);

  /* Calculate the angular acceleration. */
  /* Taken from "An improved algorithm for molecular dynamics simulation of
   * rigid molecules", Sonnenschein, Roland (1985), Eq. 5.*/
  if (p->p.rotation & ROTATION_X)
    Wd[0] = (p->f.torque[0] + p->m.omega[1] * p->m.omega[2] *
                                  (p->p.rinertia[1] - p->p.rinertia[2])) /
            p->p.rinertia[0];
  else
    Wd[0] = 0.0;
  if (p->p.rotation & ROTATION_Y)
    Wd[1] = (p->f.torque[1] + p->m.omega[2] * p->m.omega[0] *
                                  (p->p.rinertia[2] - p->p.rinertia[0])) /
            p->p.rinertia[1];
  else
    Wd[1] = 0.0;
  if (p->p.rotation & ROTATION_Z)
    Wd[2] = (p->f.torque[2] + p->m.omega[0] * p->m.omega[1] *
                                  (p->p.rinertia[0] - p->p.rinertia[1])) /
            p->p.rinertia[2];
  else
    Wd[2] = 0.0;

  auto const S1 = Qd[0] * Qd[0] + Qd[1] * Qd[1] + Qd[2] * Qd[2] + Qd[3] * Qd[3];

  /* Calculate the second derivative of the quaternion. */
  /* Taken from "An improved algorithm for molecular dynamics simulation of
   * rigid molecules", Sonnenschein, Roland (1985), Eq. 8.*/
  Qdd[0] = 0.5 * (-p->r.quat[1] * Wd[0] - p->r.quat[2] * Wd[1] -
                  p->r.quat[3] * Wd[2]) -
           p->r.quat[0] * S1;

  Qdd[1] = 0.5 * (p->r.quat[0] * Wd[0] - p->r.quat[3] * Wd[1] +
                  p->r.quat[2] * Wd[2]) -
           p->r.quat[1] * S1;

  Qdd[2] = 0.5 * (p->r.quat[3] * Wd[0] + p->r.quat[0] * Wd[1] -
                  p->r.quat[1] * Wd[2]) -
           p->r.quat[2] * S1;

  Qdd[3] = 0.5 * (-p->r.quat[2] * Wd[0] + p->r.quat[1] * Wd[1] +
                  p->r.quat[0] * Wd[2]) -
           p->r.quat[3] * S1;

  S[0] = S1;
  S[1] = Qd[0] * Qdd[0] + Qd[1] * Qdd[1] + Qd[2] * Qdd[2] + Qd[3] * Qdd[3];
  S[2] = Qdd[0] * Qdd[0] + Qdd[1] * Qdd[1] + Qdd[2] * Qdd[2] + Qdd[3] * Qdd[3];
}

/** propagate angular velocities and quaternions \todo implement for
       fixed_coord_flag */
void propagate_omega_quat_particle(Particle *p) {
  double lambda;

  double Qd[4], Qdd[4], S[3], Wd[3];
  // If rotation for the particle is disabled entirely, return early.
  if (!p->p.rotation)
    return;
#ifdef BROWNIAN_DYNAMICS
  if (thermo_switch & THERMO_BROWNIAN)
    return;
#endif // BROWNIAN_DYNAMICS

  // Clear rotational velocity for blocked rotation axes.
  if (!(p->p.rotation & ROTATION_X))
    p->m.omega[0] = 0;
  if (!(p->p.rotation & ROTATION_Y))
    p->m.omega[1] = 0;
  if (!(p->p.rotation & ROTATION_Z))
    p->m.omega[2] = 0;

  define_Qdd(p, Qd, Qdd, S, Wd);

  /* Taken from "On the numerical integration of motion for rigid polyatomics:
   * The modified quaternion approach", Omeylan, Igor (1998), Eq. 12.*/
  lambda = 1 - S[0] * time_step_squared_half -
           sqrt(1 - time_step_squared *
                        (S[0] + time_step * (S[1] + time_step_half / 2. *
                                                        (S[2] - S[0] * S[0]))));

  for (int j = 0; j < 3; j++) {
    p->m.omega[j] += time_step_half * Wd[j];
  }
  ONEPART_TRACE(if (p->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",
                            this_node, p->m.v[0], p->m.v[1], p->m.v[2]));

  p->r.quat[0] +=
      time_step * (Qd[0] + time_step_half * Qdd[0]) - lambda * p->r.quat[0];
  p->r.quat[1] +=
      time_step * (Qd[1] + time_step_half * Qdd[1]) - lambda * p->r.quat[1];
  p->r.quat[2] +=
      time_step * (Qd[2] + time_step_half * Qdd[2]) - lambda * p->r.quat[2];
  p->r.quat[3] +=
      time_step * (Qd[3] + time_step_half * Qdd[3]) - lambda * p->r.quat[3];
  // Update the director

  ONEPART_TRACE(if (p->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",
                            this_node, p->r.p[0], p->r.p[1], p->r.p[2]));
}

inline void convert_torque_to_body_frame_apply_fix_and_thermostat(Particle &p) {
  auto const t = convert_vector_space_to_body(p, p.f.torque);
  p.f.torque = Utils::Vector3d{{0, 0, 0}};

  if (thermo_switch & THERMO_LANGEVIN) {
#if defined(VIRTUAL_SITES) && defined(THERMOSTAT_IGNORE_NON_VIRTUAL)
    if (!p.p.is_virtual)
#endif
    {
      friction_thermo_langevin_rotation(&p);

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
void convert_torques_propagate_omega() {
  INTEG_TRACE(
      fprintf(stderr, "%d: convert_torques_propagate_omega:\n", this_node));

#if defined(CUDA) && defined(ENGINE)
  if ((lb_lbfluid_get_lattice_switch() == ActiveLB::GPU) &&
      swimming_particles_exist) {
    copy_v_cs_from_GPU(local_cells.particles());
  }
#endif

  for (auto &p : local_cells.particles()) {
    // Skip particle if rotation is turned off entirely for it.
    if (!p.p.rotation)
      continue;

    convert_torque_to_body_frame_apply_fix_and_thermostat(p);

#if defined(ENGINE)
    if (p.swim.swimming && lb_lbfluid_get_lattice_switch() != ActiveLB::NONE) {

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

#ifdef BROWNIAN_DYNAMICS
    if (thermo_switch & THERMO_BROWNIAN) {
      bd_drag_vel_rot(p, 0.5 * time_step);
      bd_random_walk_vel_rot(p, 0.5 * time_step);
    } else
#endif // BROWNIAN_DYNAMICS
    {
      ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
          stderr, "%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",
          this_node, p.f.f[0], p.f.f[1], p.f.f[2], p.m.v[0], p.m.v[1], p.m.v[2]));

      // Propagation of angular velocities
      p.m.omega[0] += time_step_half * p.f.torque[0] / p.p.rinertia[0];
      p.m.omega[1] += time_step_half * p.f.torque[1] / p.p.rinertia[1];
      p.m.omega[2] += time_step_half * p.f.torque[2] / p.p.rinertia[2];

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
}

/** convert the torques to the body-fixed frames before the integration loop */
void convert_initial_torques() {

  INTEG_TRACE(fprintf(stderr, "%d: convert_initial_torques:\n", this_node));
  for (auto &p : local_cells.particles()) {
    if (!p.p.rotation)
      continue;
    convert_torque_to_body_frame_apply_fix_and_thermostat(p);
  }
}
// Frame conversion routines

Utils::Vector3d convert_vector_body_to_space(const Particle &p,
                                             const Utils::Vector3d &vec) {
  Utils::Vector3d res = {0, 0, 0};
  double A[9];
  define_rotation_matrix(p, A);

  res[0] =
      A[0 + 3 * 0] * vec[0] + A[1 + 3 * 0] * vec[1] + A[2 + 3 * 0] * vec[2];
  res[1] =
      A[0 + 3 * 1] * vec[0] + A[1 + 3 * 1] * vec[1] + A[2 + 3 * 1] * vec[2];
  res[2] =
      A[0 + 3 * 2] * vec[0] + A[1 + 3 * 2] * vec[1] + A[2 + 3 * 2] * vec[2];

  return res;
}

Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                             const Utils::Vector3d &v) {
  Utils::Vector3d res = {0, 0, 0};
  double A[9];
  define_rotation_matrix(p, A);
  res[0] = A[0 + 3 * 0] * v[0] + A[0 + 3 * 1] * v[1] + A[0 + 3 * 2] * v[2];
  res[1] = A[1 + 3 * 0] * v[0] + A[1 + 3 * 1] * v[1] + A[1 + 3 * 2] * v[2];
  res[2] = A[2 + 3 * 0] * v[0] + A[2 + 3 * 1] * v[1] + A[2 + 3 * 2] * v[2];
  return res;
}

/** Fixing the per-particle per-axis rotations
 */
void rotation_fix(Particle &p, Utils::Vector3d &rot_vector) {
  // Per coordinate fixing
  if (!(p.p.rotation & ROTATION_X))
    rot_vector[0] = 0;
  if (!(p.p.rotation & ROTATION_Y))
    rot_vector[1] = 0;
  if (!(p.p.rotation & ROTATION_Z))
    rot_vector[2] = 0;
}
/** Rotate the particle p around the body-frame defined NORMALIZED axis
 * aBodyFrame by amount phi
 */
void local_rotate_particle_body(Particle &p, const Utils::Vector3d &axis_body_frame,
                           const double phi) {
  Utils::Vector3d axis = axis_body_frame;

  // Rotation turned off entirely?
  if (!p.p.rotation)
    return;

  // Per coordinate fixing
  rotation_fix(p, axis);
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

/** Rotate the particle p around the NORMALIZED axis aSpaceFrame by amount phi
 */
void local_rotate_particle(Particle &p, const Utils::Vector3d &axis_space_frame,
                           const double phi) {
  // Convert rotation axis to body-fixed frame
  Utils::Vector3d axis = convert_vector_space_to_body(p, axis_space_frame);
  local_rotate_particle_body(p, axis, phi);
}

#ifdef BROWNIAN_DYNAMICS

/** Propagate quaternions: viscous drag driven by conservative torques.*/
/*********************************************************/
/** \name bd_drag_rot */
/*********************************************************/
/**(An analogy of eq. (14.39) T. Schlick,
 * https://doi.org/10.1007/978-1-4419-6351-2 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
void bd_drag_rot(Particle &p, double dt) {
  Thermostat::GammaType local_gamma;

  if (p.p.gamma_rot >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma_rot;
  } else {
    local_gamma = langevin_gamma_rotation;
  }

  Utils::Vector3d dphi = {0.0, 0.0, 0.0};
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // only a conservative part of the torque is used here
#ifndef PARTICLE_ANISOTROPY
      dphi[j] = p.f.torque[j] * dt / (local_gamma);
#else
      dphi[j] = p.f.torque[j] * dt / (local_gamma[j]);
#endif // ROTATIONAL_INERTIA
    }
  } // j
  rotation_fix(p, dphi);
  double dphi_m = dphi.norm();
  if (dphi_m) {
    Utils::Vector3d dphi_u;
    dphi_u = dphi / dphi_m;
    local_rotate_particle_body(p, dphi_u, dphi_m);
  }
}

/** Set the terminal angular velocity driven by the conservative torques drag.*/
/*********************************************************/
/** \name bd_drag_vel_rot */
/*********************************************************/
/**(An analogy of the 1st term of the eq. (14.34) T. Schlick,
 * https://doi.org/10.1007/978-1-4419-6351-2 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
void bd_drag_vel_rot(Particle &p, double dt) {
  Thermostat::GammaType local_gamma;

  if (p.p.gamma_rot >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma_rot;
  } else {
    local_gamma = langevin_gamma_rotation;
  }

  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (p.p.ext_flag & COORD_FIXED(j)) {
      p.m.omega[j] = 0.0;
    } else
#endif
    {
      // only conservative part of the force is used here
      // NOTE: velocity is assigned here and propagated by thermal part further
      // on top of it
#ifndef PARTICLE_ANISOTROPY
      p.m.omega[j] = p.f.torque[j] / (local_gamma);
#else
      p.m.omega[j] = p.f.torque[j] / (local_gamma[j]);
#endif // ROTATIONAL_INERTIA
    }
  }
  rotation_fix(p, p.m.omega);
}

/** Propagate the quaternions: random walk part.*/
/*********************************************************/
/** \name bd_random_walk_rot */
/*********************************************************/
/**(An analogy of eq. (14.37) T. Schlick,
 * https://doi.org/10.1007/978-1-4419-6351-2 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
void bd_random_walk_rot(Particle &p, double dt) {
  extern Thermostat::GammaType brown_sigma_pos_rotation_inv;
  extern Thermostat::GammaType brown_gammatype_nan;
  // first, set defaults
  Thermostat::GammaType brown_sigma_pos_temp_inv = brown_sigma_pos_rotation_inv;

  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  auto const constexpr langevin_temp_coeff = 2.0;

  if (p.p.gamma_rot >= Thermostat::GammaType{}) {
    // Is a particle-specific temperature also specified?
    if (p.p.T >= 0.) {
      if (p.p.T > 0.0) {
        brown_sigma_pos_temp_inv =
            sqrt(p.p.gamma_rot / (langevin_temp_coeff * p.p.T));
      } else {
        brown_sigma_pos_temp_inv =
            brown_gammatype_nan; // just an indication of the infinity
      }
    } else
        // Default temperature but particle-specific gamma
        if (temperature > 0.) {
      brown_sigma_pos_temp_inv =
          sqrt(p.p.gamma_rot / (langevin_temp_coeff * temperature));
    } else {
      brown_sigma_pos_temp_inv = brown_gammatype_nan;
    }
  } // particle specific gamma
  else {
    // No particle-specific gamma, but is there particle-specific temperature
    if (p.p.T >= 0.) {
      if (p.p.T > 0.0) {
        brown_sigma_pos_temp_inv =
            sqrt(langevin_gamma_rotation / (langevin_temp_coeff * p.p.T));
      } else {
        brown_sigma_pos_temp_inv =
            brown_gammatype_nan; // just an indication of the infinity
      }
    } else {
      // Defaut values for both
      brown_sigma_pos_temp_inv = brown_sigma_pos_rotation_inv;
    }
  }
#endif /* LANGEVIN_PER_PARTICLE */

  Utils::Vector3d dphi = {0.0, 0.0, 0.0};
  Utils::Vector3d noise = v_noise_g(p.p.identity, RNGSalt::BROWNIAN);
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
#ifndef PARTICLE_ANISOTROPY
      if (brown_sigma_pos_temp_inv > 0.0) {
        dphi[j] =
            noise[j] * (1.0 / brown_sigma_pos_temp_inv) * sqrt(dt);
      } else {
        dphi[j] = 0.0;
      }
#else
      if (brown_sigma_pos_temp_inv[j] > 0.0) {
        dphi[j] = noise[j] * (1.0 / brown_sigma_pos_temp_inv[j]) *
                  sqrt(dt);
      } else {
        dphi[j] = 0.0;
      }
#endif // ROTATIONAL_INERTIA
    }
  }
  rotation_fix(p, dphi);
  // making the algorithm to be independ on an order of the rotations
  double dphi_m = dphi.norm();
  if (dphi_m) {
    Utils::Vector3d dphi_u;
    dphi_u = dphi / dphi_m;
    local_rotate_particle_body(p, dphi_u, dphi_m);
  }
}

/** Determine the angular velocities: random walk part.*/
/*********************************************************/
/** \name bd_random_walk_vel_rot */
/*********************************************************/
/**(An analogy of eq. (10.2.16) N. Pottier,
 * https://doi.org/10.1007/s10955-010-0114-6 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
void bd_random_walk_vel_rot(Particle &p, double dt) {
  extern double brown_sigma_vel_rotation;
  // first, set defaults
  double brown_sigma_vel_temp = brown_sigma_vel_rotation;

  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  auto const constexpr langevin_temp_coeff = 1.0;
  // Is a particle-specific temperature specified?
  if (p.p.T >= 0.) {
    brown_sigma_vel_temp = sqrt(langevin_temp_coeff * p.p.T);
  } else {
    brown_sigma_vel_temp = brown_sigma_vel_rotation;
  }
#endif /* LANGEVIN_PER_PARTICLE */

  Utils::Vector3d domega;
  Utils::Vector3d noise = v_noise_g(p.p.identity, RNGSalt::BROWNIAN);
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // velocity is added here. It is already initialized in the terminal drag
      // part.
      domega[j] =
          brown_sigma_vel_temp * noise[j] / sqrt(p.p.rinertia[j]);
    }
  }
  rotation_fix(p, domega);
  p.m.omega += domega;
}

#endif // BROWNIAN_DYNAMICS

#endif
