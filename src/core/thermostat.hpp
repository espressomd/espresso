/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _THERMOSTAT_H
#define _THERMOSTAT_H
/** \file thermostat.hpp

*/

#include "config.hpp"

#include <cmath>
#include "debug.hpp"
#include "particle_data.hpp"
#include "random.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "dpd.hpp"
#include "virtual_sites.hpp"
#include "Vector.hpp"

/** \name Thermostat switches*/
/************************************************************/
/*@{*/

#define THERMO_OFF        0
#define THERMO_LANGEVIN   1
#define THERMO_DPD        2
#define THERMO_NPT_ISO    4
#define THERMO_LB         8
#define THERMO_GHMC       32
/*@}*/

namespace Thermostat {
static auto noise = []() { return (d_random() - 0.5); };

#ifdef PARTICLE_ANISOTROPY
using GammaType = Vector3d;
#else
using GammaType = double;
#endif
}

/************************************************
 * exported variables
 ************************************************/

/** Switch determining which thermostat to use. This is a or'd value
    of the different possible thermostats (defines: \ref THERMO_OFF,
    \ref THERMO_LANGEVIN, \ref THERMO_DPD \ref THERMO_NPT_ISO). If it
    is zero all thermostats are switched off and the temperature is
    set to zero.  */
extern int thermo_switch;

/** temperature. */
extern double temperature;

/** Langevin friction coefficient gamma. */
extern Thermostat::GammaType langevin_gamma;
/** Langevin friction coefficient gamma. */
extern Thermostat::GammaType langevin_gamma_rotation;

/** Langevin for translations */
extern bool langevin_trans;

/** Langevin for rotations */
extern bool langevin_rotate;

/** Friction coefficient for nptiso-thermostat's inline-function
 * friction_therm0_nptiso */
extern double nptiso_gamma0;
/** Friction coefficient for nptiso-thermostat's inline-function
 * friction_thermV_nptiso */
extern double nptiso_gammav;

/** Number of NVE-MD steps in GHMC Cycle*/
extern int ghmc_nmd;
/** Phi parameter for GHMC partial momenum update step */
extern double ghmc_phi;

/************************************************
 * functions
 ************************************************/

/** initialize constants of the thermostat on
    start of integration */
void thermo_init();

/** very nasty: if we recalculate force when leaving/reentering the integrator,
    a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
    numbers are drawn twice, resulting in a different variance of the random
   force.
    This is corrected by additional heat when restarting the integrator here.
    Currently only works for the Langevin thermostat, although probably also
   others
    are affected.
*/
void thermo_heat_up();

/** pendant to \ref thermo_heat_up */
void thermo_cool_down();

#ifdef ROTATION
inline void thermo_define_rotation_matrix(Particle *p, double A[9]) {
  double q0q0 = p->r.quat[0];
  q0q0 *= q0q0;

  double q1q1 = p->r.quat[1];
  q1q1 *= q1q1;

  double q2q2 = p->r.quat[2];
  q2q2 *= q2q2;

  double q3q3 = p->r.quat[3];
  q3q3 *= q3q3;

  A[0 + 3 * 0] = q0q0 + q1q1 - q2q2 - q3q3;
  A[1 + 3 * 1] = q0q0 - q1q1 + q2q2 - q3q3;
  A[2 + 3 * 2] = q0q0 - q1q1 - q2q2 + q3q3;

  A[0 + 3 * 1] =
      2 * (p->r.quat[1] * p->r.quat[2] + p->r.quat[0] * p->r.quat[3]);
  A[0 + 3 * 2] =
      2 * (p->r.quat[1] * p->r.quat[3] - p->r.quat[0] * p->r.quat[2]);
  A[1 + 3 * 0] =
      2 * (p->r.quat[1] * p->r.quat[2] - p->r.quat[0] * p->r.quat[3]);

  A[1 + 3 * 2] =
      2 * (p->r.quat[2] * p->r.quat[3] + p->r.quat[0] * p->r.quat[1]);
  A[2 + 3 * 0] =
      2 * (p->r.quat[1] * p->r.quat[3] + p->r.quat[0] * p->r.quat[2]);
  A[2 + 3 * 1] =
      2 * (p->r.quat[2] * p->r.quat[3] - p->r.quat[0] * p->r.quat[1]);
}

inline void thermo_convert_forces_body_to_space(Particle *p, double *force) {
  double A[9];
  thermo_define_rotation_matrix(p, A);

  force[0] = A[0 + 3 * 0] * p->f.f[0] + A[1 + 3 * 0] * p->f.f[1] +
             A[2 + 3 * 0] * p->f.f[2];
  force[1] = A[0 + 3 * 1] * p->f.f[0] + A[1 + 3 * 1] * p->f.f[1] +
             A[2 + 3 * 1] * p->f.f[2];
  force[2] = A[0 + 3 * 2] * p->f.f[0] + A[1 + 3 * 2] * p->f.f[1] +
             A[2 + 3 * 2] * p->f.f[2];
}

inline void thermo_convert_vel_space_to_body(Particle *p, double *vel_space,
                                             double *vel_body) {
  double A[9];
  thermo_define_rotation_matrix(p, A);

  vel_body[0] = A[0 + 3 * 0] * vel_space[0] + A[0 + 3 * 1] * vel_space[1] +
                A[0 + 3 * 2] * vel_space[2];
  vel_body[1] = A[1 + 3 * 0] * vel_space[0] + A[1 + 3 * 1] * vel_space[1] +
                A[1 + 3 * 2] * vel_space[2];
  vel_body[2] = A[2 + 3 * 0] * vel_space[0] + A[2 + 3 * 1] * vel_space[1] +
                A[2 + 3 * 2] * vel_space[2];
}
#endif // ROTATION

/** locally defined funcion to find Vx. In case of LEES_EDWARDS, that is
   relative to the LE shear frame
    @param i      coordinate index
    @param vel    velocity vector
    @param pos    position vector
    @return       adjusted (or not) i^th velocity coordinate */
inline double le_frameV(int i, double *vel, double *pos) {
#ifdef LEES_EDWARDS

  if (i == 0) {
    double relY = pos[1] * box_l_i[1] - 0.5;
    return (vel[0] - relY * lees_edwards_rate);
  }

#endif

  return vel[i];
}

#ifdef NPT
/** add velocity-dependend noise and friction for NpT-sims to the particle's
   velocity
    @param dt_vj  j-component of the velocity scaled by time_step dt
    @return       j-component of the noise added to the velocity, also scaled by
   dt (contained in prefactors) */
inline double friction_therm0_nptiso(double dt_vj) {
  extern double nptiso_pref1, nptiso_pref2;
  if (thermo_switch & THERMO_NPT_ISO)
    return (nptiso_pref1 * dt_vj + nptiso_pref2 * Thermostat::noise());

  return 0.0;
}

/** add p_diff-dependend noise and friction for NpT-sims to \ref
 * nptiso_struct::p_diff */
inline double friction_thermV_nptiso(double p_diff) {
  extern double nptiso_pref3, nptiso_pref4;
  if (thermo_switch & THERMO_NPT_ISO)
    return (nptiso_pref3 * p_diff + nptiso_pref4 * Thermostat::noise());
  return 0.0;
}
#endif

/** overwrite the forces of a particle with
    the friction term, i.e. \f$ F_i= -\gamma v_i + \xi_i\f$.
*/
inline void friction_thermo_langevin(Particle *p) {
  extern Thermostat::GammaType langevin_pref1, langevin_pref2;
  Thermostat::GammaType langevin_pref1_temp, langevin_pref2_temp;

#ifdef MULTI_TIMESTEP
  extern double langevin_pref1_small;
#ifndef LANGEVIN_PER_PARTICLE
  extern double langevin_pref2_small;
#endif /* LANGEVIN_PER_PARTICLE */
#endif /* MULTI_TIMESTEP */

  int j;
  double switch_trans = 1.0;
  if (langevin_trans == false) {
    switch_trans = 0.0;
  }

// Virtual sites related decision making
#ifdef VIRTUAL_SITES
#ifndef VIRTUAL_SITES_THERMOSTAT
  // In this case, virtual sites are NOT thermostated
  if (p->p.isVirtual) {
    for (j = 0; j < 3; j++)
      p->f.f[j] = 0;

    return;
  }
#endif /* VIRTUAL_SITES_THERMOSTAT */
#ifdef THERMOSTAT_IGNORE_NON_VIRTUAL
  // In this case NON-virtual particles are NOT thermostated
  if (!p->p.isVirtual) {
    for (j = 0; j < 3; j++)
      p->f.f[j] = 0;

    return;
  }
#endif /* THERMOSTAT_IGNORE_NON_VIRTUAL */
#endif /* VIRTUAL_SITES */

  // Get velocity effective in the thermostatting
  double velocity[3];
  for (int i = 0; i < 3; i++) {
    // Particle velocity
    velocity[i] = p->m.v[i];
#ifdef ENGINE
    // In case of the engine feature, the velocity is relaxed
    // towards a swimming velocity oriented parallel to the
    // particles director
    velocity[i] -= (p->swim.v_swim * time_step) * p->r.quatu[i];
#endif

    // Local effective velocity for leeds-edwards boundary conditions
    velocity[i] = le_frameV(i, velocity, p->r.p);
  } // for

// Determine prefactors for the friction and the noise term

// first, set defaults
  langevin_pref1_temp = langevin_pref1;
  langevin_pref2_temp = langevin_pref2;

// Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  auto const constexpr langevin_temp_coeff = 24.0;

  if (p->p.gamma >= Thermostat::GammaType{}) {
    langevin_pref1_temp = -p->p.gamma / time_step;
    // Is a particle-specific temperature also specified?
    if (p->p.T >= 0.)
      langevin_pref2_temp =
          sqrt(langevin_temp_coeff * p->p.T * p->p.gamma / time_step);
    else
      // Default temperature but particle-specific gamma
      langevin_pref2_temp =
          sqrt(langevin_temp_coeff * temperature * p->p.gamma / time_step);

  } // particle specific gamma
  else {
    langevin_pref1_temp = -langevin_gamma / time_step;
    // No particle-specific gamma, but is there particle-specific temperature
    if (p->p.T >= 0.)
      langevin_pref2_temp =
          sqrt(langevin_temp_coeff * p->p.T * langevin_gamma / time_step);
    else
      // Defaut values for both
      langevin_pref2_temp = langevin_pref2;
  }

#endif /* LANGEVIN_PER_PARTICLE */

// Multi-timestep handling
// This has to be last, as it may set the prefactors to 0.
#ifdef MULTI_TIMESTEP
  if (smaller_time_step > 0.) {
    langevin_pref1_temp *= time_step / smaller_time_step;
    if (p->p.smaller_timestep == 1 && current_time_step_is_small == 1)
      langevin_pref2_temp *= sqrt(time_step / smaller_time_step);
    else if (p->p.smaller_timestep != current_time_step_is_small) {
      langevin_pref1_temp = 0.;
      langevin_pref2_temp = 0.;
    }
  }
#endif /* MULTI_TIMESTEP */

#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check
  auto aniso_flag = (langevin_pref1_temp[0] != langevin_pref1_temp[1]) ||
               (langevin_pref1_temp[1] != langevin_pref1_temp[2]) ||
               (langevin_pref2_temp[0] != langevin_pref2_temp[1]) ||
               (langevin_pref2_temp[1] != langevin_pref2_temp[2]);
  double velocity_body[3] = {0.0, 0.0, 0.0};
  if (aniso_flag) {
     thermo_convert_vel_space_to_body(p, velocity, velocity_body);
  }
#endif

  // Do the actual thermostatting
  for (j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    // If individual coordinates are fixed, set force to 0.
    if ((p->p.ext_flag & COORD_FIXED(j)))
      p->f.f[j] = 0;
    else
#endif
    {
// Apply the force
#ifndef PARTICLE_ANISOTROPY
      p->f.f[j] = langevin_pref1_temp * velocity[j] +
                  switch_trans * langevin_pref2_temp * Thermostat::noise();
#else
      // In case of anisotropic particle: body-fixed reference frame. Otherwise:
      // lab-fixed reference frame.
      if (aniso_flag)
        p->f.f[j] = langevin_pref1_temp[j] * velocity_body[j] +
                    switch_trans * langevin_pref2_temp[j] * Thermostat::noise();
      else
        p->f.f[j] = langevin_pref1_temp[j] * velocity[j] +
                    switch_trans * langevin_pref2_temp[j] * Thermostat::noise();
#endif
    }
  } // END LOOP OVER ALL COMPONENTS

#ifdef PARTICLE_ANISOTROPY
  if (aniso_flag) {
    double particle_force[3] = {0.0, 0.0, 0.0};

    thermo_convert_forces_body_to_space(p, particle_force);
    for (j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p->p.ext_flag & COORD_FIXED(j)))
#endif
      {
        p->f.f[j] = particle_force[j];
      }
    }
  }
#endif // PARTICLE_ANISOTROPY

  // printf("%d: %e %e %e %e %e %e\n",p->p.identity,
  // p->f.f[0],p->f.f[1],p->f.f[2], p->m.v[0],p->m.v[1],p->m.v[2]);
  ONEPART_TRACE(if (p->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",
                            this_node, p->f.f[0], p->f.f[1], p->f.f[2]));
  THERMO_TRACE(fprintf(stderr, "%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",
                       this_node, p->p.identity, p->f.f[0], p->f.f[1],
                       p->f.f[2]));
}

#ifdef ROTATION
/** set the particle torques to the friction term, i.e. \f$\tau_i=-\gamma w_i +
   \xi_i\f$.
    The same friction coefficient \f$\gamma\f$ is used as that for translation.
*/
inline void friction_thermo_langevin_rotation(Particle *p) {
  extern Thermostat::GammaType langevin_pref2_rotation;
  Thermostat::GammaType langevin_pref1_temp, langevin_pref2_temp;

  double switch_rotate = 1.0;
  if (langevin_rotate == false) {
    switch_rotate = 0.0;
  }

  langevin_pref1_temp = langevin_gamma_rotation;
  langevin_pref2_temp = langevin_pref2_rotation;

// Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  // If a particle-specific gamma is given
  auto const constexpr langevin_temp_coeff = 24.0;

  if (p->p.gamma_rot >= Thermostat::GammaType{}) {
    langevin_pref1_temp = p->p.gamma_rot;
    // Is a particle-specific temperature also specified?
    if (p->p.T >= 0.)
      langevin_pref2_temp =
          sqrt(langevin_temp_coeff * p->p.T * p->p.gamma_rot / time_step);
    else
      // Default temperature but particle-specific gamma
      langevin_pref2_temp =
          sqrt(langevin_temp_coeff * temperature * p->p.gamma_rot / time_step);

  } // particle specific gamma
  else {
    langevin_pref1_temp = langevin_gamma_rotation;
    // No particle-specific gamma, but is there particle-specific temperature
    if (p->p.T >= 0.)
      langevin_pref2_temp = sqrt(langevin_temp_coeff * p->p.T *
                                 langevin_gamma_rotation / time_step);
    else
      // Default values for both
      langevin_pref2_temp = langevin_pref2_rotation;
  }
#endif /* LANGEVIN_PER_PARTICLE */

  // Rotational degrees of virtual sites are thermostatted,
  // so no switching here

  // Here the thermostats happens
  for (int j = 0; j < 3; j++) {
#ifdef PARTICLE_ANISOTROPY
    p->f.torque[j] =
        -langevin_pref1_temp[j] * p->m.omega[j] +
        switch_rotate * langevin_pref2_temp[j] * Thermostat::noise();
#else
    p->f.torque[j] = -langevin_pref1_temp * p->m.omega[j] +
                     switch_rotate * langevin_pref2_temp * Thermostat::noise();
#endif
  }

  ONEPART_TRACE(if (p->p.identity == check_id)
                    fprintf(stderr, "%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",
                            this_node, p->f.f[0], p->f.f[1], p->f.f[2]));
  THERMO_TRACE(fprintf(stderr, "%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",
                       this_node, p->p.identity, p->f.f[0], p->f.f[1],
                       p->f.f[2]));
}

#endif // ROTATION
#endif
