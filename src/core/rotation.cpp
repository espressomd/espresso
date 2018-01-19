/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file rotation.cpp  Molecular dynamics integrator for rotational motion.
 *
 *  A velocity Verlet <a
 * HREF="http://ciks.cbt.nist.gov/~garbocz/dpd1/dpd.html">algorithm</a>
 *  using quaternions is implemented to tackle rotational motion. A random
 * torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature
 * all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of
 * freedom if ROTATION
 *  flag is set in \ref config.hpp "config.h".
*/

#include "rotation.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "forces.hpp"
#include "ghosts.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "p3m.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"
#include "utils.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

/****************************************************
 *                     DEFINES
 ***************************************************/
/**************** local variables  *******************/

#ifdef ROTATION

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** define first and second time derivatives of a quaternion, as well
    as the angular acceleration. */
static void define_Qdd(Particle *p, double Qd[4], double Qdd[4], double S[3],
                       double Wd[3]);

/*@}*/

/** convert quaternions to the director */
/** Convert director to quaternions */
int convert_quatu_to_quat(double d[3], double quat[4]) {
  double d_xy, dm;
  double theta2, phi2;

  // Calculate magnitude of the given vector
  dm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);

  // The vector needs to be != 0 to be converted into a quaternion
  if (dm < ROUND_ERROR_PREC) {
    return 1;
  } else {
    // Calculate angles
    d_xy = sqrt(d[0] * d[0] + d[1] * d[1]);
    // If dipole points along z axis:
    if (d_xy == 0) {
      // We need to distinguish between (0,0,d_z) and (0,0,d_z)
      if (d[2] > 0)
        theta2 = 0;
      else
        theta2 = PI / 2.;
      phi2 = 0;
    } else {
      // Here, we take care of all other directions
      // Here we suppose that theta2 = 0.5*theta and phi2 = 0.5*(phi - PI/2),
      // where theta and phi - angles are in spherical coordinates
      theta2 = 0.5 * acos(d[2] / dm);
      if (d[1] < 0)
        phi2 = -0.5 * acos(d[0] / d_xy) - PI * 0.25;
      else
        phi2 = 0.5 * acos(d[0] / d_xy) - PI * 0.25;
    }

    // Calculate the quaternion from the angles
    quat[0] = cos(theta2) * cos(phi2);
    quat[1] = -sin(theta2) * cos(phi2);
    quat[2] = -sin(theta2) * sin(phi2);
    quat[3] = cos(theta2) * sin(phi2);
  }
  return 0;
}

/** Here we use quaternions to calculate the rotation matrix which
    will be used then to transform torques from the laboratory to
    the body-fixed frames */  
void define_rotation_matrix(Particle const *p, double A[9])
{
  double q0q0 =p->r.quat[0];
  q0q0 *=q0q0;

  double q1q1 =p->r.quat[1];
  q1q1 *=q1q1;

  double q2q2 =p->r.quat[2];
  q2q2 *=q2q2;

  double q3q3 =p->r.quat[3];
  q3q3 *=q3q3;

  A[0 + 3*0] = q0q0 + q1q1 - q2q2 - q3q3;
  A[1 + 3*1] = q0q0 - q1q1 + q2q2 - q3q3;
  A[2 + 3*2] = q0q0 - q1q1 - q2q2 + q3q3;

  A[0 + 3*1] = 2*(p->r.quat[1]*p->r.quat[2] + p->r.quat[0]*p->r.quat[3]);
  A[0 + 3*2] = 2*(p->r.quat[1]*p->r.quat[3] - p->r.quat[0]*p->r.quat[2]);
  A[1 + 3*0] = 2*(p->r.quat[1]*p->r.quat[2] - p->r.quat[0]*p->r.quat[3]);

  A[1 + 3*2] = 2*(p->r.quat[2]*p->r.quat[3] + p->r.quat[0]*p->r.quat[1]);
  A[2 + 3*0] = 2*(p->r.quat[1]*p->r.quat[3] + p->r.quat[0]*p->r.quat[2]);
  A[2 + 3*1] = 2*(p->r.quat[2]*p->r.quat[3] - p->r.quat[0]*p->r.quat[1]);
}

/** calculate the second derivative of the quaternion of a given particle
    as well as Wd vector which is the angular acceleration of this particle */
void define_Qdd(Particle *p, double Qd[4], double Qdd[4], double S[3],
                double Wd[3]) {
  /* calculate the first derivative of the quaternion */
  Qd[0] = 0.5 * (-p->r.quat[1] * p->m.omega[0] - p->r.quat[2] * p->m.omega[1] -
                 p->r.quat[3] * p->m.omega[2]);

  Qd[1] = 0.5 * (p->r.quat[0] * p->m.omega[0] - p->r.quat[3] * p->m.omega[1] +
                 p->r.quat[2] * p->m.omega[2]);

  Qd[2] = 0.5 * (p->r.quat[3] * p->m.omega[0] + p->r.quat[0] * p->m.omega[1] -
                 p->r.quat[1] * p->m.omega[2]);

  Qd[3] = 0.5 * (-p->r.quat[2] * p->m.omega[0] + p->r.quat[1] * p->m.omega[1] +
                 p->r.quat[0] * p->m.omega[2]);

  /* calculate the second derivative of the quaternion */  
  Wd[0] =  (p->f.torque[0] + p->m.omega[1]*p->m.omega[2]*(p->p.rinertia[1]-p->p.rinertia[2]))/p->p.rinertia[0];
  Wd[1] =  (p->f.torque[1] + p->m.omega[2]*p->m.omega[0]*(p->p.rinertia[2]-p->p.rinertia[0]))/p->p.rinertia[1];
  Wd[2] =  (p->f.torque[2] + p->m.omega[0]*p->m.omega[1]*(p->p.rinertia[0]-p->p.rinertia[1]))/p->p.rinertia[2];

  auto const S1 = Qd[0] * Qd[0] + Qd[1] * Qd[1] + Qd[2] * Qd[2] + Qd[3] * Qd[3];

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

  // Clear rotational velocity for blocked rotation axes.
  if (! (p->p.rotation & ROTATION_X))
    p->m.omega[0]=0;
  if (! (p->p.rotation & ROTATION_Y))
    p->m.omega[1]=0;
  if (! (p->p.rotation & ROTATION_Z))
    p->m.omega[2]=0;


  
  define_Qdd(p, Qd, Qdd, S, Wd);

  lambda = 1 - S[0] * time_step_squared_half -
           sqrt(1 -
                time_step_squared *
                    (S[0] +
                     time_step *
                         (S[1] + time_step_half / 2. * (S[2] - S[0] * S[0]))));

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
  convert_quat_to_quatu(p->r.quat, p->r.quatu);
#ifdef DIPOLES
  // When dipoles are enabled, update dipole moment
  convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif

  ONEPART_TRACE(if (p->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",
                            this_node, p->r.p[0], p->r.p[1], p->r.p[2]));
}

/** convert the torques to the body-fixed frames and propagate angular
 * velocities */
void convert_torques_propagate_omega() {
  double tx, ty, tz;
  double omega_0[3] = {0.0, 0.0, 0.0};

  INTEG_TRACE(
      fprintf(stderr, "%d: convert_torques_propagate_omega:\n", this_node));

#if defined(LB_GPU) && defined(ENGINE)
  if (lattice_switch & LATTICE_LB_GPU) {
    copy_v_cs_from_GPU(local_cells.particles());
  }
#endif

  for (auto &p : local_cells.particles()) {
    // Skip particle if rotation is turned off entirely for it.
    if (!p.p.rotation)
      continue;
    
    double A[9];
    define_rotation_matrix(&p, A);

    tx = A[0 + 3 * 0] * p.f.torque[0] + A[0 + 3 * 1] * p.f.torque[1] +
         A[0 + 3 * 2] * p.f.torque[2];
    ty = A[1 + 3 * 0] * p.f.torque[0] + A[1 + 3 * 1] * p.f.torque[1] +
         A[1 + 3 * 2] * p.f.torque[2];
    tz = A[2 + 3 * 0] * p.f.torque[0] + A[2 + 3 * 1] * p.f.torque[1] +
         A[2 + 3 * 2] * p.f.torque[2];

    if (thermo_switch & THERMO_LANGEVIN) {
#if defined(VIRTUAL_SITES) && defined(THERMOSTAT_IGNORE_NON_VIRTUAL)
      if (!p.p.isVirtual)
#endif
      {
        friction_thermo_langevin_rotation(&p);

        p.f.torque[0] += tx;
        p.f.torque[1] += ty;
        p.f.torque[2] += tz;
      }
    } else {
      p.f.torque[0] = tx;
      p.f.torque[1] = ty;
      p.f.torque[2] = tz;
    }

    if (!(p.p.rotation & ROTATION_X))
      p.f.torque[0] = 0;

    if (!(p.p.rotation & ROTATION_Y))
      p.f.torque[1] = 0;

    if (!(p.p.rotation & ROTATION_Z))
      p.f.torque[2] = 0;


#if defined(ENGINE) && (defined(LB) || defined(LB_GPU))
    double omega_swim[3] = {0, 0, 0};
    double omega_swim_body[3] = {0, 0, 0};
    if (p.swim.swimming && lattice_switch != 0) {
      double dip[3];
      double diff[3];
      double cross[3];
      double l_diff, l_cross;

      dip[0] = p.swim.dipole_length * p.r.quatu[0];
      dip[1] = p.swim.dipole_length * p.r.quatu[1];
      dip[2] = p.swim.dipole_length * p.r.quatu[2];

      diff[0] = (p.swim.v_center[0] - p.swim.v_source[0]);
      diff[1] = (p.swim.v_center[1] - p.swim.v_source[1]);
      diff[2] = (p.swim.v_center[2] - p.swim.v_source[2]);

      cross[0] = diff[1] * dip[2] - diff[2] * dip[1];
      cross[1] = diff[0] * dip[2] - diff[2] * dip[0];
      cross[2] = diff[0] * dip[1] - diff[1] * dip[0];

      l_diff = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      l_cross =
          sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

      if (l_cross > 0 && p.swim.dipole_length > 0) {
        omega_swim[0] = l_diff * cross[0] / (l_cross * p.swim.dipole_length);
        omega_swim[1] = l_diff * cross[1] / (l_cross * p.swim.dipole_length);
        omega_swim[2] = l_diff * cross[2] / (l_cross * p.swim.dipole_length);

        omega_swim_body[0] = A[0 + 3 * 0] * omega_swim[0] +
                             A[0 + 3 * 1] * omega_swim[1] +
                             A[0 + 3 * 2] * omega_swim[2];
        omega_swim_body[1] = A[1 + 3 * 0] * omega_swim[0] +
                             A[1 + 3 * 1] * omega_swim[1] +
                             A[1 + 3 * 2] * omega_swim[2];
        omega_swim_body[2] = A[2 + 3 * 0] * omega_swim[0] +
                             A[2 + 3 * 1] * omega_swim[1] +
                             A[2 + 3 * 2] * omega_swim[2];

        p.f.torque[0] +=
            p.swim.rotational_friction * (omega_swim_body[0] - p.m.omega[0]);
        p.f.torque[1] +=
            p.swim.rotational_friction * (omega_swim_body[1] - p.m.omega[1]);
        p.f.torque[2] +=
            p.swim.rotational_friction * (omega_swim_body[2] - p.m.omega[2]);
      }
    }
#endif

    ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
        stderr, "%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",
        this_node, p.f.f[0], p.f.f[1], p.f.f[2], p.m.v[0], p.m.v[1], p.m.v[2]));

    p.m.omega[0] += time_step_half * p.f.torque[0] / p.p.rinertia[0];
    p.m.omega[1] += time_step_half * p.f.torque[1] / p.p.rinertia[1];
    p.m.omega[2] += time_step_half * p.f.torque[2] / p.p.rinertia[2];

    // zeroth estimate of omega
    for (int j = 0; j < 3; j++)
      omega_0[j] = p.m.omega[j];

    /* if the tensor of inertia is isotropic, the following refinement is not
       needed.
       Otherwise repeat this loop 2-3 times depending on the required accuracy
       */
    for (int times = 0; times <= 5; times++) {
      double Wd[3];

      Wd[0] =
          (p.m.omega[1] * p.m.omega[2] * (p.p.rinertia[1] - p.p.rinertia[2])) /
          p.p.rinertia[0];
      Wd[1] =
          (p.m.omega[2] * p.m.omega[0] * (p.p.rinertia[2] - p.p.rinertia[0])) /
          p.p.rinertia[1];
      Wd[2] =
          (p.m.omega[0] * p.m.omega[1] * (p.p.rinertia[0] - p.p.rinertia[1])) /
          p.p.rinertia[2];

      p.m.omega[0] = omega_0[0] + time_step_half * Wd[0];
      p.m.omega[1] = omega_0[1] + time_step_half * Wd[1];
      p.m.omega[2] = omega_0[2] + time_step_half * Wd[2];
    }

    ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
        stderr, "%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n", this_node, p.m.v[0],
        p.m.v[1], p.m.v[2]));
  }
}

/** convert the torques to the body-fixed frames before the integration loop */
void convert_initial_torques() {
  double tx, ty, tz;

  INTEG_TRACE(fprintf(stderr, "%d: convert_initial_torques:\n", this_node));
  for (auto &p : local_cells.particles()) {
    if (!p.p.rotation)
      continue;
    double A[9];
    define_rotation_matrix(&p, A);

    tx = A[0 + 3 * 0] * p.f.torque[0] + A[0 + 3 * 1] * p.f.torque[1] +
         A[0 + 3 * 2] * p.f.torque[2];
    ty = A[1 + 3 * 0] * p.f.torque[0] + A[1 + 3 * 1] * p.f.torque[1] +
         A[1 + 3 * 2] * p.f.torque[2];
    tz = A[2 + 3 * 0] * p.f.torque[0] + A[2 + 3 * 1] * p.f.torque[1] +
         A[2 + 3 * 2] * p.f.torque[2];

    if (thermo_switch & THERMO_LANGEVIN) {

      friction_thermo_langevin_rotation(&p);
      p.f.torque[0] += tx;
      p.f.torque[1] += ty;
      p.f.torque[2] += tz;
    } else {
      p.f.torque[0] = tx;
      p.f.torque[1] = ty;
      p.f.torque[2] = tz;
    }

    if (!(p.p.rotation & ROTATION_X))
      p.f.torque[0] = 0;

    if (!(p.p.rotation & ROTATION_Y))
      p.f.torque[1] = 0;

    if (!(p.p.rotation & ROTATION_Z))
      p.f.torque[2] = 0;


    ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
        stderr, "%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",
        this_node, p.f.f[0], p.f.f[1], p.f.f[2], p.m.v[0], p.m.v[1], p.m.v[2]));
  }
}

/** convert from the body-fixed frames to space-fixed coordinates */

void convert_omega_body_to_space(const Particle *p, double *omega) {
  double A[9];
  define_rotation_matrix(p, A);

  omega[0] = A[0 + 3 * 0] * p->m.omega[0] + A[1 + 3 * 0] * p->m.omega[1] +
             A[2 + 3 * 0] * p->m.omega[2];
  omega[1] = A[0 + 3 * 1] * p->m.omega[0] + A[1 + 3 * 1] * p->m.omega[1] +
             A[2 + 3 * 1] * p->m.omega[2];
  omega[2] = A[0 + 3 * 2] * p->m.omega[0] + A[1 + 3 * 2] * p->m.omega[1] +
             A[2 + 3 * 2] * p->m.omega[2];
}

Vector3d convert_vector_body_to_space(const Particle& p, const Vector3d& vec) {
  Vector3d res={0,0,0};
  double A[9];
  define_rotation_matrix(&p, A);

  res[0] = A[0 + 3 * 0] * vec[0] + A[1 + 3 * 0] * vec[1] +
             A[2 + 3 * 0] * vec[2];
  res[1] = A[0 + 3 * 1] * vec[0] + A[1 + 3 * 1] * vec[1] +
             A[2 + 3 * 1] * vec[2];
  res[2] = A[0 + 3 * 2] * vec[0] + A[1 + 3 * 2] * vec[1] +
             A[2 + 3 * 2] * vec[2];
  
  return res;
}


void convert_torques_body_to_space(const Particle *p, double *torque) {
  double A[9];
  define_rotation_matrix(p, A);

  torque[0] = A[0 + 3 * 0] * p->f.torque[0] + A[1 + 3 * 0] * p->f.torque[1] +
              A[2 + 3 * 0] * p->f.torque[2];
  torque[1] = A[0 + 3 * 1] * p->f.torque[0] + A[1 + 3 * 1] * p->f.torque[1] +
              A[2 + 3 * 1] * p->f.torque[2];
  torque[2] = A[0 + 3 * 2] * p->f.torque[0] + A[1 + 3 * 2] * p->f.torque[1] +
              A[2 + 3 * 2] * p->f.torque[2];
}

void convert_vel_space_to_body(const Particle *p, double *vel_body) {
  double A[9];
  define_rotation_matrix(p, A);

  vel_body[0] = A[0 + 3 * 0] * p->m.v[0] + A[0 + 3 * 1] * p->m.v[1] +
                A[0 + 3 * 2] * p->m.v[2];
  vel_body[1] = A[1 + 3 * 0] * p->m.v[0] + A[1 + 3 * 1] * p->m.v[1] +
                A[1 + 3 * 2] * p->m.v[2];
  vel_body[2] = A[2 + 3 * 0] * p->m.v[0] + A[2 + 3 * 1] * p->m.v[1] +
                A[2 + 3 * 2] * p->m.v[2];
}

void convert_vec_space_to_body(Particle *p, double *v, double *res) {
  double A[9];
  define_rotation_matrix(p, A);

  res[0] = A[0 + 3 * 0] * v[0] + A[0 + 3 * 1] * v[1] + A[0 + 3 * 2] * v[2];
  res[1] = A[1 + 3 * 0] * v[0] + A[1 + 3 * 1] * v[1] + A[1 + 3 * 2] * v[2];
  res[2] = A[2 + 3 * 0] * v[0] + A[2 + 3 * 1] * v[1] + A[2 + 3 * 2] * v[2];
}


/** Rotate the particle p around the NORMALIZED axis aSpaceFrame by amount phi
 */
void rotate_particle(Particle *p, double *aSpaceFrame, double phi) {
  // Convert rotation axis to body-fixed frame
  double a[3];
  convert_vec_space_to_body(p, aSpaceFrame, a);

  //  printf("%g %g %g - ",a[0],a[1],a[2]);
  // Rotation turned off entirely?
  if (!p->p.rotation)
    return;

  // Per coordinate fixing
  if (!(p->p.rotation & ROTATION_X))
    a[0] = 0;
  if (!(p->p.rotation & ROTATION_Y))
    a[1] = 0;
  if (!(p->p.rotation & ROTATION_Z))
    a[2] = 0;
  // Re-normalize rotation axis
  double l = sqrt(sqrlen(a));
  // Check, if the rotation axis is nonzero
  if (l < 1E-10)
    return;

  for (int i = 0; i < 3; i++)
    a[i] /= l;

  double q[4];
  q[0] = cos(phi / 2);
  double tmp = sin(phi / 2);
  q[1] = tmp * a[0];
  q[2] = tmp * a[1];
  q[3] = tmp * a[2];

  // Normalize
  normalize_quaternion(q);

  // Rotate the particle
  double qn[4]; // Resulting quaternion
  multiply_quaternions(p->r.quat, q, qn);
  for (int k = 0; k < 4; k++)
    p->r.quat[k] = qn[k];
  convert_quat_to_quatu(p->r.quat, p->r.quatu);
#ifdef DIPOLES
  // When dipoles are enabled, update dipole moment
  convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
}

#endif
