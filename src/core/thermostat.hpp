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
#ifndef CORE_THERMOSTAT_HPP
#define CORE_THERMOSTAT_HPP
/** \file

*/

#include "config.hpp"

#include "debug.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "random.hpp"
#include "rotation.hpp"

#include "utils/Vector.hpp"

#include <cmath>

/** \name Thermostat switches*/
/************************************************************/
/*@{*/

#define THERMO_OFF 0
#define THERMO_LANGEVIN 1
#define THERMO_DPD 2
#define THERMO_NPT_ISO 4
#define THERMO_LB 8
/*@}*/

namespace Thermostat {
static auto noise = []() { return (d_random() - 0.5); };

#ifdef PARTICLE_ANISOTROPY
using GammaType = Vector3d;
#else
using GammaType = double;
#endif
} // namespace Thermostat

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

/** True if the thermostat should act on virtual particles. */
extern bool thermo_virtual;

/** Langevin friction coefficient gamma. */
extern Thermostat::GammaType langevin_gamma;
/** Langevin friction coefficient gamma. */
extern Thermostat::GammaType langevin_gamma_rotation;

/** Friction coefficient for nptiso-thermostat's inline-function
 * friction_therm0_nptiso */
extern double nptiso_gamma0;
/** Friction coefficient for nptiso-thermostat's inline-function
 * friction_thermV_nptiso */
extern double nptiso_gammav;

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

#ifdef NPT
/** add velocity-dependent noise and friction for NpT-sims to the particle's
   velocity
    @param vj     j-component of the velocity
    @return       j-component of the noise added to the velocity, also scaled by
   dt (contained in prefactors) */
inline double friction_therm0_nptiso(double vj) {
  extern double nptiso_pref1, nptiso_pref2;
  if (thermo_switch & THERMO_NPT_ISO) {
    if (nptiso_pref2 > 0.0) {
      return (nptiso_pref1 * vj + nptiso_pref2 * Thermostat::noise());
    } else {
      return nptiso_pref1 * vj;
    }
  }
  return 0.0;
}

/** add p_diff-dependent noise and friction for NpT-sims to \ref
 * nptiso_struct::p_diff */
inline double friction_thermV_nptiso(double p_diff) {
  extern double nptiso_pref3, nptiso_pref4;
  if (thermo_switch & THERMO_NPT_ISO) {
    if (nptiso_pref4 > 0.0) {
      return (nptiso_pref3 * p_diff + nptiso_pref4 * Thermostat::noise());
    } else {
      return nptiso_pref3 * p_diff;
    }
  }
  return 0.0;
}
#endif

/** overwrite the forces of a particle with
    the friction term, i.e. \f$ F_i= -\gamma v_i + \xi_i\f$.
*/
inline void friction_thermo_langevin(Particle *p) {
  extern Thermostat::GammaType langevin_pref1, langevin_pref2;
  Thermostat::GammaType langevin_pref1_temp, langevin_pref2_temp;

  if (p->p.is_virtual && !thermo_virtual) {
    for (int j = 0; j < 3; j++)
      p->f.f[j] = 0;

    return;
  }

  // Get velocity effective in the thermostatting
  Vector3d velocity = p->m.v;
#ifdef ENGINE
  if (p->swim.v_swim != 0) {
    // In case of the engine feature, the velocity is relaxed
    // towards a swimming velocity oriented parallel to the
    // particles director
    velocity -= p->swim.v_swim * p->r.calc_director();
  }
#endif

  // Determine prefactors for the friction and the noise term

  // first, set defaults
  langevin_pref1_temp = langevin_pref1;
  langevin_pref2_temp = langevin_pref2;

// Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  auto const constexpr langevin_temp_coeff = 24.0;

  if (p->p.gamma >= Thermostat::GammaType{}) {
    langevin_pref1_temp = -p->p.gamma;
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
    langevin_pref1_temp = -langevin_gamma;
    // No particle-specific gamma, but is there particle-specific temperature
    if (p->p.T >= 0.)
      langevin_pref2_temp =
          sqrt(langevin_temp_coeff * p->p.T * langevin_gamma / time_step);
    else
      // Default values for both
      langevin_pref2_temp = langevin_pref2;
  }

#endif /* LANGEVIN_PER_PARTICLE */

#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check
  auto aniso_flag = (langevin_pref1_temp[0] != langevin_pref1_temp[1]) ||
                    (langevin_pref1_temp[1] != langevin_pref1_temp[2]) ||
                    (langevin_pref2_temp[0] != langevin_pref2_temp[1]) ||
                    (langevin_pref2_temp[1] != langevin_pref2_temp[2]);
  auto const velocity_body =
      (aniso_flag) ? convert_vector_space_to_body(*p, velocity) : Vector3d{};
#endif

  // Do the actual thermostatting
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    // If individual coordinates are fixed, set force to 0.
    if ((p->p.ext_flag & COORD_FIXED(j)))
      p->f.f[j] = 0;
    else
#endif
    {
// Apply the force
#ifndef PARTICLE_ANISOTROPY
      if (langevin_pref2_temp > 0.0) {
        p->f.f[j] = langevin_pref1_temp * velocity[j] +
                    langevin_pref2_temp * Thermostat::noise();
      } else {
        p->f.f[j] = langevin_pref1_temp * velocity[j];
      }
#else
      // In case of anisotropic particle: body-fixed reference frame. Otherwise:
      // lab-fixed reference frame.
      if (aniso_flag) {
        if (langevin_pref2_temp[j] > 0.0) {
          p->f.f[j] = langevin_pref1_temp[j] * velocity_body[j] +
                      langevin_pref2_temp[j] * Thermostat::noise();
        } else {
          p->f.f[j] = langevin_pref1_temp[j] * velocity_body[j];
        }
      } else {
        if (langevin_pref2_temp[j] > 0.0) {
          p->f.f[j] = langevin_pref1_temp[j] * velocity[j] +
                      langevin_pref2_temp[j] * Thermostat::noise();
        } else {
          p->f.f[j] = langevin_pref1_temp[j] * velocity[j];
        }
      }
#endif
    }
  } // END LOOP OVER ALL COMPONENTS

#ifdef PARTICLE_ANISOTROPY
  if (aniso_flag) {
    auto const particle_force = convert_vector_body_to_space(*p, p->f.f);

    for (int j = 0; j < 3; j++) {
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
    if (langevin_pref2_temp[j] > 0.0) {
      p->f.torque[j] = -langevin_pref1_temp[j] * p->m.omega[j] +
                       langevin_pref2_temp[j] * Thermostat::noise();
    } else {
      p->f.torque[j] = -langevin_pref1_temp[j] * p->m.omega[j];
    }
#else
    if (langevin_pref2_temp > 0.0) {
      p->f.torque[j] = -langevin_pref1_temp * p->m.omega[j] +
                       langevin_pref2_temp * Thermostat::noise();
    } else {
      p->f.torque[j] = -langevin_pref1_temp * p->m.omega[j];
    }
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
