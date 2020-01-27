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
  /** Recalculate prefactors.
   *  Needs to be called every time the parameters are changed.
   */
  void recalc_prefactors() {
    pref_friction = -gamma;
    pref_noise = sigma(temperature, time_step, gamma);
    // If gamma_rotation is not set explicitly, use the translational one.
    if (gamma_rotation < GammaType{}) {
      gamma_rotation = gamma;
    }
    pref_noise_rotation = sigma(temperature, time_step, gamma_rotation);
  }
  /** Calculate the noise standard deviation. */
  static GammaType sigma(double kT, double time_step, GammaType const &gamma) {
    constexpr auto const temp_coeff = 24.0;
    return sqrt((temp_coeff * kT / time_step) * gamma);
  }
  /** @name Parameters */
  /*@{*/
  /** Translational friction coefficient @f$ \gamma_{\text{trans}} @f$. */
  GammaType gamma = sentinel(GammaType{});
  /** Rotational friction coefficient @f$ \gamma_{\text{rot}} @f$. */
  GammaType gamma_rotation = sentinel(GammaType{});
  /*@}*/
  /** @name Prefactors */
  /*@{*/
  /** Prefactor for the friction. */
  GammaType pref_friction;
  /** Prefactor for the translational velocity noise. */
  GammaType pref_noise;
  /** Prefactor for the angular velocity noise. */
  GammaType pref_noise_rotation;
  /*@}*/
  /** RNG counter, used for both translation and rotation. */
  std::unique_ptr<Utils::Counter<uint64_t>> rng_counter;
};

/** %Thermostat for Brownian dynamics.
 *  Default particle mass is assumed to be unitary in these global parameters.
 */
struct BrownianThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Recalculate prefactors.
   *  Needs to be called every time the parameters are changed.
   */
  void recalc_prefactors() {
    /** The heat velocity dispersion corresponds to the Gaussian noise only,
     *  which is only valid for the BD. Just a square root of kT, see (10.2.17)
     *  and comments in 2 paragraphs afterwards, @cite Pottier2010.
     */
    sigma_vel = sigma(temperature);
    /** The random walk position dispersion is defined by the second eq. (14.38)
     *  of @cite Schlick2010. Its time interval factor will be added in the
     *  Brownian Dynamics functions. Its square root is the standard deviation.
     */
    sigma_pos = sigma(temperature, gamma);
#ifdef ROTATION
    /** Note: the BD thermostat assigns the brownian viscous parameters as well.
     *  They correspond to the friction tensor Z from the eq. (14.31) of
     *  @cite Schlick2010.
     */
    // If gamma_rotation is not set explicitly, use the translational one.
    if (gamma_rotation < GammaType{}) {
      gamma_rotation = gamma;
    }
    sigma_vel_rotation = sigma(temperature);
    sigma_pos_rotation = sigma(temperature, gamma_rotation);
#endif // ROTATION
  }
  /** Calculate the noise standard deviation. */
  static GammaType sigma(double kT, GammaType const &gamma) {
    constexpr auto const temp_coeff = 2.0;
    return sqrt(Utils::hadamard_division(temp_coeff * kT, gamma));
  }
  /** Calculate the noise standard deviation. */
  static double sigma(double kT) {
    constexpr auto const temp_coeff = 1.0;
    return sqrt(temp_coeff * kT);
  }
  /** @name Parameters */
  /*@{*/
  /** Translational friction coefficient @f$ \gamma_{\text{trans}} @f$. */
  GammaType gamma = sentinel(GammaType{});
  /** Rotational friction coefficient @f$ \gamma_{\text{rot}} @f$. */
  GammaType gamma_rotation = sentinel(GammaType{});
  /*@}*/
  /** @name Prefactors */
  /*@{*/
  /** Translational noise standard deviation.
   *  Stores @f$ \sqrt{2D_{\text{trans}}} @f$ with
   *  @f$ D_{\text{trans}} = k_B T/\gamma_{\text{trans}} @f$
   *  the translational diffusion coefficient.
   */
  GammaType sigma_pos = sentinel(GammaType{});
  /** Rotational noise standard deviation.
   *  Stores @f$ \sqrt{2D_{\text{rot}}} @f$ with
   *  @f$ D_{\text{rot}} = k_B T/\gamma_{\text{rot}} @f$
   *  the rotational diffusion coefficient.
   */
  GammaType sigma_pos_rotation = sentinel(GammaType{});
  /** Translational velocity noise standard deviation.
   *  Stores @f$ \sqrt{k_B T} @f$.
   */
  double sigma_vel = 0;
  /** Angular velocity noise standard deviation.
   *  Stores @f$ \sqrt{k_B T} @f$.
   */
  double sigma_vel_rotation = 0;
  /*@}*/
  /** RNG counter, used for both translation and rotation. */
  std::unique_ptr<Utils::Counter<uint64_t>> rng_counter;
};

/** %Thermostat for isotropic NPT dynamics. */
struct IsotropicNptThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
#ifdef NPT
  /** Recalculate prefactors.
   *  Needs to be called every time the parameters are changed.
   */
  void recalc_prefactors(double piston) {
    assert(piston > 0.0);
    pref1 = -gamma0 * 0.5 * time_step;
    pref2 = sqrt(12.0 * temperature * gamma0 * time_step);
    pref3 = -gammav * (1.0 / piston) * 0.5 * time_step;
    pref4 = sqrt(12.0 * temperature * gammav * time_step);
  }
#endif
  /** @name Parameters */
  /*@{*/
  /** Friction coefficient @f$ \gamma_0 @f$ */
  double gamma0;
  /** Friction coefficient @f$ \gamma_V @f$ */
  double gammav;
  /*@}*/
#ifdef NPT
  /** @name Prefactors */
  /*@{*/
  double pref1;
  double pref2;
  double pref3;
  double pref4;
  /*@}*/
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
