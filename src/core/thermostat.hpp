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
#define THERMO_BROWNIAN 16
/*@}*/

namespace Thermostat {

static auto noise = []() { return (d_random() - 0.5); };

#ifdef PARTICLE_ANISOTROPY
using GammaType = Utils::Vector3d;
#else
using GammaType = double;
#endif
} // namespace Thermostat

namespace {
/** @name Integrators parameters sentinels.
 *  These functions return the sentinel value for the Langevin/Brownian
 *  parameters, indicating that they have not been set yet.
 */
/*@{*/
constexpr double sentinel(double) { return -1.0; }
constexpr Utils::Vector3d sentinel(Utils::Vector3d) {
  return {-1.0, -1.0, -1.0};
}
constexpr double set_nan(double) { return NAN; }
constexpr Utils::Vector3d set_nan(Utils::Vector3d) { return {NAN, NAN, NAN}; }
/*@}*/
} // namespace

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

/************************************************
 * parameter structs
 ************************************************/

/** %Thermostat for Langevin dynamics. */
struct LangevinThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Translational friction coefficient @f$ \gamma_{\text{trans}} @f$. */
  GammaType gamma = sentinel(GammaType{});
  /** Rotational friction coefficient @f$ \gamma_{\text{rot}} @f$. */
  GammaType gamma_rotation = sentinel(GammaType{});
  /** Prefactor for the friction. */
  GammaType pref_friction;
  /** Prefactor for the translational velocity noise. */
  GammaType pref_noise;
  /** Prefactor for the angular velocity noise. */
  GammaType pref_noise_rotation;
  /** RNG counter, used for both translation and rotation. */
  std::unique_ptr<Utils::Counter<uint64_t>> rng_counter;
};

/** %Thermostat for Brownian dynamics. */
struct BrownianThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Translational friction coefficient @f$ \gamma_{\text{trans}} @f$. */
  GammaType gamma = sentinel(GammaType{});
  /** Rotational friction coefficient @f$ \gamma_{\text{rot}} @f$. */
  GammaType gamma_rotation = sentinel(GammaType{});
  /** Inverse of the translational noise standard deviation.
   *  Stores @f$ \left(\sqrt{2D_{\text{trans}}}\right)^{-1} @f$ with
   *  @f$ D_{\text{trans}} = k_B T/\gamma_{\text{trans}} @f$
   *  the translational diffusion coefficient
   */
  GammaType sigma_pos_inv = sentinel(GammaType{});
  /** Inverse of the rotational noise standard deviation.
   *  Stores @f$ \left(\sqrt{2D_{\text{rot}}}\right)^{-1} @f$ with
   *  @f$ D_{\text{rot}} = k_B T/\gamma_{\text{rot}} @f$
   *  the rotational diffusion coefficient
   */
  GammaType sigma_pos_rotation_inv = sentinel(GammaType{});
  /** Sentinel value for divisions by zero. */
  GammaType const gammatype_nan = set_nan(GammaType{});
  /** Translational velocity noise standard deviation. */
  double sigma_vel = 0;
  /** Angular velocity noise standard deviation. */
  double sigma_vel_rotation = 0;
  /** RNG counter, used for both translation and rotation. */
  std::unique_ptr<Utils::Counter<uint64_t>> rng_counter;
};

/** %Thermostat for isotropic NPT dynamics. */
struct IsotropicNptThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Friction coefficient @f$ \gamma_0 @f$ */
  double gamma0;
  /** Friction coefficient @f$ \gamma_V @f$ */
  double gammav;
#ifdef NPT
  double pref1;
  double pref2;
  double pref3;
  double pref4;
#endif
};

/************************************************
 * functions
 ************************************************/

/** Only require seed if rng is not initialized. */
bool langevin_is_seed_required();

/** Only require seed if rng is not initialized. */
bool brownian_is_seed_required();

/** @name philox functionality: increment, get/set */
/*@{*/
void langevin_rng_counter_increment();
void langevin_set_rng_state(uint64_t counter);
uint64_t langevin_get_rng_state();
void brownian_rng_counter_increment();
void brownian_set_rng_state(uint64_t counter);
uint64_t brownian_get_rng_state();
/*@}*/

/** Initialize constants of the thermostat at the start of integration */
void thermo_init();

#ifdef NPT
/** Add velocity-dependent noise and friction for NpT-sims to the particle's
 *  velocity
 *  @param npt_iso Parameters
 *  @param vj     j-component of the velocity
 *  @return       j-component of the noise added to the velocity, also scaled by
 *                dt (contained in prefactors)
 */
inline double friction_therm0_nptiso(IsotropicNptThermostat const &npt_iso,
                                     double vj) {
  if (thermo_switch & THERMO_NPT_ISO) {
    if (npt_iso.pref2 > 0.0) {
      return npt_iso.pref1 * vj + npt_iso.pref2 * Thermostat::noise();
    }
    return npt_iso.pref1 * vj;
  }
  return 0.0;
}

/** Add p_diff-dependent noise and friction for NpT-sims to \ref
 *  nptiso_struct::p_diff
 */
inline double friction_thermV_nptiso(IsotropicNptThermostat const &npt_iso,
                                     double p_diff) {
  if (thermo_switch & THERMO_NPT_ISO) {
    if (npt_iso.pref4 > 0.0) {
      return npt_iso.pref3 * p_diff + npt_iso.pref4 * Thermostat::noise();
    }
    return npt_iso.pref3 * p_diff;
  }
  return 0.0;
}
#endif

#endif
