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

#include <boost/optional.hpp>
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

struct BaseThermostat {
public:
  /** Initialize or re-initialize the RNG counter with a seed. */
  void rng_initialize(uint64_t const seed) {
    rng_counter = Utils::Counter<uint64_t>(seed);
  }
  /** Increment the RNG counter */
  void rng_increment() {
    if (!rng_counter) {
      throw "The RNG counter is not initialized";
    }
    rng_counter.get().increment();
  }
  /** Get current value of the RNG */
  uint64_t rng_get() const {
    if (!rng_counter) {
      throw "The RNG counter is not initialized";
    }
    return rng_counter.get().value();
  }
  /** Is the RNG counter initialized */
  bool rng_is_initialized() const { return static_cast<bool>(rng_counter); }

private:
  /** RNG counter. */
  boost::optional<Utils::Counter<uint64_t>> rng_counter;
};

/** %Thermostat for Langevin dynamics. */
struct LangevinThermostat : public BaseThermostat {
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
  /** Calculate the noise prefactor.
   *  Evaluates the quantity @f$ \sqrt{2 k_B T \gamma / dt} / \sigma_\eta @f$
   *  with @f$ \sigma_\eta @f$ the standard deviation of the random uniform
   *  process @f$ \eta(t) @f$.
   */
  static GammaType sigma(double kT, double time_step, GammaType const &gamma) {
    // random uniform noise has variance 1/12
    constexpr auto const temp_coeff = 2.0 * 12.0;
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
  /** Prefactor for the friction.
   *  Stores @f$ \gamma_{\text{trans}} @f$.
   */
  GammaType pref_friction;
  /** Prefactor for the translational velocity noise.
   *  Stores @f$ \sqrt{2 k_B T \gamma_{\text{trans}} / dt} / \sigma_\eta @f$.
   */
  GammaType pref_noise;
  /** Prefactor for the angular velocity noise.
   *  Stores @f$ \sqrt{2 k_B T \gamma_{\text{rot}} / dt} / \sigma_\eta @f$.
   */
  GammaType pref_noise_rotation;
  /*@}*/
};

/** %Thermostat for Brownian dynamics.
 *  Default particle mass is assumed to be unitary in these global parameters.
 */
struct BrownianThermostat : public BaseThermostat {
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
  /** Calculate the noise prefactor.
   *  Evaluates the quantity @f$ \sqrt{2 k_B T / \gamma} / \sigma_\eta @f$
   *  with @f$ \sigma_\eta @f$ the standard deviation of the random gaussian
   *  process @f$ \eta(t) @f$.
   */
  static GammaType sigma(double kT, GammaType const &gamma) {
    constexpr auto const temp_coeff = 2.0;
    return sqrt(Utils::hadamard_division(temp_coeff * kT, gamma));
  }
  /** Calculate the noise prefactor.
   *  Evaluates the quantity @f$ \sqrt{k_B T} / \sigma_\eta @f$
   *  with @f$ \sigma_\eta @f$ the standard deviation of the random gaussian
   *  process @f$ \eta(t) @f$.
   */
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
};

/** %Thermostat for isotropic NPT dynamics. */
struct IsotropicNptThermostat : public BaseThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Recalculate prefactors.
   *  Needs to be called every time the parameters are changed.
   */
  void recalc_prefactors(double piston) {
#ifdef NPT
    assert(piston > 0.0);
    auto const half_time_step = time_step / 2.0;
    pref_rescale_0 = -gamma0 * half_time_step;
    pref_noise_0 = sigma(temperature, gamma0);
    pref_rescale_V = -gammav * half_time_step / piston;
    pref_noise_V = sigma(temperature, gammav);
#endif
  }
  /** Calculate the noise prefactor.
   *  Evaluates the quantity @f$ \sqrt{2 k_B T \gamma dt / 2} / \sigma_\eta @f$
   *  with @f$ \sigma_\eta @f$ the standard deviation of the random uniform
   *  process @f$ \eta(t) @f$.
   */
  static double sigma(double kT, double gamma) {
    // random uniform noise has variance 1/12; the temperature
    // coefficient of 2 is canceled out by the half time step
    constexpr auto const temp_coeff = 12.0;
    return sqrt(temp_coeff * temperature * gamma * time_step);
  }
  /** @name Parameters */
  /*@{*/
  /** Friction coefficient of the particles @f$ \gamma^0 @f$ */
  double gamma0;
  /** Friction coefficient for the box @f$ \gamma^V @f$ */
  double gammav;
  /*@}*/
#ifdef NPT
  /** @name Prefactors */
  /*@{*/
  /** %Particle velocity rescaling at half the time step.
   *  Stores @f$ \gamma^{0}\cdot\frac{dt}{2} @f$.
   */
  double pref_rescale_0;
  /** %Particle velocity rescaling noise standard deviation.
   *  Stores @f$ \sqrt{k_B T \gamma^{0} dt} / \sigma_\eta @f$.
   */
  double pref_noise_0;
  /** Volume rescaling at half the time step.
   *  Stores @f$ \frac{\gamma^{V}}{Q}\cdot\frac{dt}{2} @f$.
   */
  double pref_rescale_V;
  /** Volume rescaling noise standard deviation.
   *  Stores @f$ \sqrt{k_B T \gamma^{V} dt} / \sigma_\eta @f$.
   */
  double pref_noise_V;
  /*@}*/
#endif
};

/** %Thermostat for thermalized bonds. */
struct ThermalizedBondThermostat : public BaseThermostat {};

#ifdef DPD
/** %Thermostat for dissipative particle dynamics. */
struct DPDThermostat : public BaseThermostat {};
#endif

/************************************************
 * functions
 ************************************************/

/**
 * @brief Register a thermostat public interface
 *
 * @param thermostat        The thermostat name
 */
#define NEW_THERMOSTAT(thermostat)                                             \
  bool thermostat##_is_seed_required();                                        \
  void thermostat##_rng_counter_increment();                                   \
  void thermostat##_set_rng_state(uint64_t counter);                           \
  uint64_t thermostat##_get_rng_state();

NEW_THERMOSTAT(langevin)
NEW_THERMOSTAT(brownian)
NEW_THERMOSTAT(npt_iso)
NEW_THERMOSTAT(thermalized_bond)
#ifdef DPD
NEW_THERMOSTAT(dpd)
#endif

/** Initialize constants of the thermostat at the start of integration */
void thermo_init();

#ifdef NPT
/** Add velocity-dependent noise and friction for NpT-sims to the particle's
 *  velocity
 *  @tparam step       Which half time step to integrate (1 or 2)
 *  @param npt_iso     Parameters
 *  @param vel         particle velocity
 *  @param p_identity  particle identity
 *  @return noise added to the velocity, already rescaled by
 *          dt/2 (contained in prefactors)
 */
template <size_t step>
inline Utils::Vector3d
friction_therm0_nptiso(IsotropicNptThermostat const &npt_iso,
                       Utils::Vector3d const &vel, int p_identity) {
  static_assert(step == 1 or step == 2, "NPT only has 2 integration steps");
  constexpr auto const salt =
      (step == 1) ? RNGSalt::NPTISO0_HALF_STEP1 : RNGSalt::NPTISO0_HALF_STEP2;
  if (thermo_switch & THERMO_NPT_ISO) {
    if (npt_iso.pref_noise_0 > 0.0) {
      return npt_iso.pref_rescale_0 * vel +
             npt_iso.pref_noise_0 *
                 Random::noise_uniform<salt>(npt_iso.rng_get(), p_identity);
    }
    return npt_iso.pref_rescale_0 * vel;
  }
  return {};
}

/** Add p_diff-dependent noise and friction for NpT-sims to \ref
 *  nptiso_struct::p_diff
 */
inline double friction_thermV_nptiso(IsotropicNptThermostat const &npt_iso,
                                     double p_diff) {
  if (thermo_switch & THERMO_NPT_ISO) {
    if (npt_iso.pref_noise_V > 0.0) {
      return npt_iso.pref_rescale_V * p_diff +
             npt_iso.pref_noise_V * Random::noise_uniform<RNGSalt::NPTISOV, 1>(
                                        npt_iso.rng_get(), 0);
    }
    return npt_iso.pref_rescale_V * p_diff;
  }
  return 0.0;
}
#endif // NPT

#endif
