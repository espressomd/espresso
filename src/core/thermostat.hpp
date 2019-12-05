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
#ifndef CORE_THERMOSTAT_HPP
#define CORE_THERMOSTAT_HPP
/** \file
 *  Implementation in \ref thermostat.cpp.
 */

#include "config.hpp"

#include "Particle.hpp"
#include "integrate.hpp"
#include "random.hpp"
#include "rotation.hpp"

#include <utils/Counter.hpp>
#include <utils/Vector.hpp>
#include <utils/math/rotation_matrix.hpp>

#include <cmath>
#include <tuple>

/** \name Thermostat switches */
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
using GammaType = Utils::Vector3d;
#else
using GammaType = double;
#endif
} // namespace Thermostat

/************************************************
 * exported variables
 ************************************************/

/** Switch determining which thermostat(s) to use. This is a or'd value
 *  of the different possible thermostats (defines: \ref THERMO_OFF,
 *  \ref THERMO_LANGEVIN, \ref THERMO_DPD \ref THERMO_NPT_ISO). If it
 *  is zero all thermostats are switched off and the temperature is
 *  set to zero.
 */
extern int thermo_switch;

/** Temperature of the thermostat. */
extern double temperature;

/** True if the thermostat should act on virtual particles. */
extern bool thermo_virtual;

/** Langevin friction coefficient gamma for translation. */
extern Thermostat::GammaType langevin_gamma;
/** Langevin friction coefficient gamma for rotation. */
extern Thermostat::GammaType langevin_gamma_rotation;

/** Friction coefficient for nptiso-thermostat's function
 *  @ref friction_therm0_nptiso
 */
extern double nptiso_gamma0;
/** Friction coefficient for nptiso-thermostat's function
 *  @ref friction_thermV_nptiso
 */
extern double nptiso_gammav;

/** Langevin RNG counter, used for both translation and rotation. */
extern std::unique_ptr<Utils::Counter<uint64_t>> langevin_rng_counter;

/************************************************
 * functions
 ************************************************/

/** Only require seed if rng is not initialized. */
bool langevin_is_seed_required();

/** @name philox functionality: increment, get/set */
/*@{*/
void langevin_rng_counter_increment();
void langevin_set_rng_state(uint64_t counter);
uint64_t langevin_get_rng_state();
/*@}*/

/** Initialize constants of the thermostat at the start of integration */
void thermo_init();

#ifdef NPT
/** Add velocity-dependent noise and friction for NpT-sims to the particle's
 *  velocity
 *  @param vj     j-component of the velocity
 *  @return       j-component of the noise added to the velocity, also scaled by
 *                dt (contained in prefactors)
 */
inline double friction_therm0_nptiso(double vj) {
  extern double nptiso_pref1, nptiso_pref2;
  if (thermo_switch & THERMO_NPT_ISO) {
    if (nptiso_pref2 > 0.0) {
      return (nptiso_pref1 * vj + nptiso_pref2 * Thermostat::noise());
    }
    return nptiso_pref1 * vj;
  }
  return 0.0;
}

/** Add p_diff-dependent noise and friction for NpT-sims to \ref
 *  nptiso_struct::p_diff
 */
inline double friction_thermV_nptiso(double p_diff) {
  extern double nptiso_pref3, nptiso_pref4;
  if (thermo_switch & THERMO_NPT_ISO) {
    if (nptiso_pref4 > 0.0) {
      return (nptiso_pref3 * p_diff + nptiso_pref4 * Thermostat::noise());
    }
    return nptiso_pref3 * p_diff;
  }
  return 0.0;
}
#endif

/** Langevin thermostat for particle translational velocities.
 *  Collects the particle velocity (different for ENGINE, PARTICLE_ANISOTROPY).
 *  Collects the langevin parameters kT, gamma (different for
 *  LANGEVIN_PER_PARTICLE). Applies the noise and friction term.
 */
inline Utils::Vector3d friction_thermo_langevin(Particle const &p) {
  // Early exit for virtual particles without thermostat
  if (p.p.is_virtual && !thermo_virtual) {
    return {};
  }

  // Determine prefactors for the friction (pref1) and the noise (pref2) term
  extern Thermostat::GammaType langevin_pref1;
  extern Thermostat::GammaType langevin_pref2;
  // first, set defaults
  Thermostat::GammaType langevin_pref_friction_buf = langevin_pref1;
  Thermostat::GammaType langevin_pref_noise_buf = langevin_pref2;
  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  if (p.p.gamma >= Thermostat::GammaType{} or p.p.T >= 0.) {
    auto const constexpr langevin_temp_coeff = 24.0;
    auto const kT = p.p.T >= 0. ? p.p.T : temperature;
    auto const gamma =
        p.p.gamma >= Thermostat::GammaType{} ? p.p.gamma : langevin_gamma;
    langevin_pref_friction_buf = -gamma;
    langevin_pref_noise_buf =
        sqrt(langevin_temp_coeff * kT * gamma / time_step);
  }
#endif /* LANGEVIN_PER_PARTICLE */

  // Get effective velocity in the thermostatting
#ifdef ENGINE
  auto const &velocity = (p.p.swim.v_swim != 0)
                             ? p.m.v - p.p.swim.v_swim * p.r.calc_director()
                             : p.m.v;
#else
  auto const &velocity = p.m.v;
#endif
#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check
  auto const aniso_flag =
      (langevin_pref_friction_buf[0] != langevin_pref_friction_buf[1]) ||
      (langevin_pref_friction_buf[1] != langevin_pref_friction_buf[2]);

  // In case of anisotropic particle: body-fixed reference frame. Otherwise:
  // lab-fixed reference frame.
  auto const friction_op =
      aniso_flag
          ? convert_body_to_space(p, diag_matrix(langevin_pref_friction_buf))
          : diag_matrix(langevin_pref_friction_buf);
  auto const noise_op = diag_matrix(langevin_pref_noise_buf);
#else
  auto const &friction_op = langevin_pref_friction_buf;
  auto const &noise_op = langevin_pref_noise_buf;
#endif // PARTICLE_ANISOTROPY
  using Random::v_noise;

  // Do the actual (isotropic) thermostatting
  return friction_op * velocity +
         noise_op * v_noise<RNGSalt::LANGEVIN>(langevin_rng_counter->value(),
                                               p.p.identity);
}

#ifdef ROTATION
/** Langevin thermostat for particle angular velocities.
 *  Collects the particle velocity (different for PARTICLE_ANISOTROPY).
 *  Collects the langevin parameters kT, gamma_rot (different for
 *  LANGEVIN_PER_PARTICLE). Applies the noise and friction term.
 */
inline Utils::Vector3d friction_thermo_langevin_rotation(const Particle &p) {
  extern Thermostat::GammaType langevin_pref2_rotation;

  auto langevin_pref_friction_buf = -langevin_gamma_rotation;
  auto langevin_pref_noise_buf = langevin_pref2_rotation;

  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  if (p.p.gamma_rot >= Thermostat::GammaType{} or p.p.T >= 0.) {
    auto const constexpr langevin_temp_coeff = 24.0;
    auto const kT = p.p.T >= 0. ? p.p.T : temperature;
    auto const gamma = p.p.gamma_rot >= Thermostat::GammaType{}
                           ? p.p.gamma_rot
                           : langevin_gamma_rotation;
    langevin_pref_friction_buf = -gamma;
    langevin_pref_noise_buf =
        sqrt(langevin_temp_coeff * kT * gamma / time_step);
  }
#endif /* LANGEVIN_PER_PARTICLE */

  using Random::v_noise;

  // Here the thermostats happens
  auto const noise = v_noise<RNGSalt::LANGEVIN_ROT>(
      langevin_rng_counter->value(), p.p.identity);
#ifdef PARTICLE_ANISOTROPY
  return hadamard_product(langevin_pref_friction_buf, p.m.omega) +
         hadamard_product(langevin_pref_noise_buf, noise);
#else
  return langevin_pref_friction_buf * p.m.omega +
         langevin_pref_noise_buf * noise;
#endif
}

#endif // ROTATION
#endif
