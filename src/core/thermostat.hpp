/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "Particle.hpp"
#include "rotation.hpp"

#include "config/config.hpp"

#include <utils/Counter.hpp>
#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>

/** \name Thermostat switches */
/**@{*/
#define THERMO_OFF 0
#define THERMO_LANGEVIN 1
#define THERMO_DPD 2
#define THERMO_NPT_ISO 4
#define THERMO_LB 8
#define THERMO_BROWNIAN 16
#define THERMO_SD 32
/**@}*/

namespace Thermostat {
#ifdef PARTICLE_ANISOTROPY
using GammaType = Utils::Vector3d;
#else
using GammaType = double;
#endif
} // namespace Thermostat

namespace {
/**
 * @brief Value for unset friction coefficient.
 * Sentinel value for the Langevin/Brownian parameters,
 * indicating that they have not been set yet.
 */
#ifdef PARTICLE_ANISOTROPY
constexpr Thermostat::GammaType gamma_sentinel{{-1.0, -1.0, -1.0}};
#else
constexpr Thermostat::GammaType gamma_sentinel{-1.0};
#endif
/**
 * @brief Value for a null friction coefficient.
 */
#ifdef PARTICLE_ANISOTROPY
constexpr Thermostat::GammaType gamma_null{{0.0, 0.0, 0.0}};
#else
constexpr Thermostat::GammaType gamma_null{0.0};
#endif
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

#ifdef THERMOSTAT_PER_PARTICLE
inline auto const &
handle_particle_gamma(Thermostat::GammaType const &particle_gamma,
                      Thermostat::GammaType const &default_gamma) {
  return particle_gamma >= gamma_null ? particle_gamma : default_gamma;
}
#endif

inline auto
handle_particle_anisotropy(Particle const &p,
                           Thermostat::GammaType const &gamma_body) {
#ifdef PARTICLE_ANISOTROPY
  auto const aniso_flag =
      (gamma_body[0] != gamma_body[1]) || (gamma_body[1] != gamma_body[2]);
  const Utils::Matrix<double, 3, 3> gamma_matrix =
      boost::qvm::diag_mat(gamma_body);
  auto const gamma_space =
      aniso_flag ? convert_body_to_space(p, gamma_matrix) : gamma_matrix;
  return gamma_space;
#else
  return gamma_body;
#endif
}

/************************************************
 * parameter structs
 ************************************************/

struct BaseThermostat {
public:
  /** Initialize or re-initialize the RNG counter with a seed. */
  void rng_initialize(uint32_t const seed) { m_rng_seed = seed; }
  /** Increment the RNG counter */
  void rng_increment() { m_rng_counter.increment(); }
  /** Get current value of the RNG */
  uint64_t rng_counter() const { return m_rng_counter.value(); }
  void set_rng_counter(uint64_t value) {
    m_rng_counter = Utils::Counter<uint64_t>(0u, value);
  }
  /** Is the RNG seed required */
  bool is_seed_required() const { return !m_rng_seed; }
  uint32_t rng_seed() const { return m_rng_seed.value(); }

private:
  /** RNG counter. */
  Utils::Counter<uint64_t> m_rng_counter;
  /** RNG seed */
  boost::optional<uint32_t> m_rng_seed;
};

/** %Thermostat for Langevin dynamics. */
struct LangevinThermostat : public BaseThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Recalculate prefactors.
   *  Needs to be called every time the parameters are changed.
   */
  void recalc_prefactors(double kT, double time_step) {
    pref_friction = -gamma;
    pref_noise = sigma(kT, time_step, gamma);
    // If gamma_rotation is not set explicitly, use the translational one.
    if (gamma_rotation < GammaType{}) {
      gamma_rotation = gamma;
    }
    pref_noise_rotation = sigma(kT, time_step, gamma_rotation);
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
  /**@{*/
  /** Translational friction coefficient @f$ \gamma_{\text{trans}} @f$. */
  GammaType gamma = gamma_sentinel;
  /** Rotational friction coefficient @f$ \gamma_{\text{rot}} @f$. */
  GammaType gamma_rotation = gamma_sentinel;
  /**@}*/
  /** @name Prefactors */
  /**@{*/
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
  /**@}*/
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
  void recalc_prefactors(double kT) {
    /** The heat velocity dispersion corresponds to the Gaussian noise only,
     *  which is only valid for the BD. Just a square root of kT, see (10.2.17)
     *  and comments in 2 paragraphs afterwards, @cite pottier10a.
     */
    sigma_vel = sigma(kT);
    /** The random walk position dispersion is defined by the second eq. (14.38)
     *  of @cite schlick10a. Its time interval factor will be added in the
     *  Brownian Dynamics functions. Its square root is the standard deviation.
     */
    sigma_pos = sigma(kT, gamma);
#ifdef ROTATION
    /** Note: the BD thermostat assigns the brownian viscous parameters as well.
     *  They correspond to the friction tensor Z from the eq. (14.31) of
     *  @cite schlick10a.
     */
    // If gamma_rotation is not set explicitly, use the translational one.
    if (gamma_rotation < GammaType{}) {
      gamma_rotation = gamma;
    }
    sigma_vel_rotation = sigma(kT);
    sigma_pos_rotation = sigma(kT, gamma_rotation);
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
  /**@{*/
  /** Translational friction coefficient @f$ \gamma_{\text{trans}} @f$. */
  GammaType gamma = gamma_sentinel;
  /** Rotational friction coefficient @f$ \gamma_{\text{rot}} @f$. */
  GammaType gamma_rotation = gamma_sentinel;
  /**@}*/
  /** @name Prefactors */
  /**@{*/
  /** Translational noise standard deviation.
   *  Stores @f$ \sqrt{2D_{\text{trans}}} @f$ with
   *  @f$ D_{\text{trans}} = k_B T/\gamma_{\text{trans}} @f$
   *  the translational diffusion coefficient.
   */
  GammaType sigma_pos = gamma_sentinel;
  /** Rotational noise standard deviation.
   *  Stores @f$ \sqrt{2D_{\text{rot}}} @f$ with
   *  @f$ D_{\text{rot}} = k_B T/\gamma_{\text{rot}} @f$
   *  the rotational diffusion coefficient.
   */
  GammaType sigma_pos_rotation = gamma_sentinel;
  /** Translational velocity noise standard deviation.
   *  Stores @f$ \sqrt{k_B T} @f$.
   */
  double sigma_vel = 0;
  /** Angular velocity noise standard deviation.
   *  Stores @f$ \sqrt{k_B T} @f$.
   */
  double sigma_vel_rotation = 0;
  /**@}*/
};

#ifdef NPT
/** %Thermostat for isotropic NPT dynamics. */
struct IsotropicNptThermostat : public BaseThermostat {
private:
  using GammaType = Thermostat::GammaType;

public:
  /** Recalculate prefactors.
   *  Needs to be called every time the parameters are changed.
   */
  void recalc_prefactors(double kT, double piston, double time_step) {
    assert(piston > 0.0);
    auto const half_time_step = time_step / 2.0;
    pref_rescale_0 = -gamma0 * half_time_step;
    pref_noise_0 = sigma(kT, gamma0, time_step);
    pref_rescale_V = -gammav * half_time_step / piston;
    pref_noise_V = sigma(kT, gammav, time_step);
  }
  /** Calculate the noise prefactor.
   *  Evaluates the quantity @f$ \sqrt{2 k_B T \gamma dt / 2} / \sigma_\eta @f$
   *  with @f$ \sigma_\eta @f$ the standard deviation of the random uniform
   *  process @f$ \eta(t) @f$.
   */
  static double sigma(double kT, double gamma, double time_step) {
    // random uniform noise has variance 1/12; the temperature
    // coefficient of 2 is canceled out by the half time step
    constexpr auto const temp_coeff = 12.0;
    return sqrt(temp_coeff * kT * gamma * time_step);
  }
  /** @name Parameters */
  /**@{*/
  /** Friction coefficient of the particles @f$ \gamma^0 @f$ */
  double gamma0;
  /** Friction coefficient for the box @f$ \gamma^V @f$ */
  double gammav;
  /**@}*/
  /** @name Prefactors */
  /**@{*/
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
  /**@}*/
};
#endif

/** %Thermostat for thermalized bonds. */
struct ThermalizedBondThermostat : public BaseThermostat {};

#ifdef DPD
/** %Thermostat for dissipative particle dynamics. */
struct DPDThermostat : public BaseThermostat {};
#endif

#ifdef STOKESIAN_DYNAMICS
/** %Thermostat for Stokesian dynamics. */
struct StokesianThermostat : public BaseThermostat {
  StokesianThermostat() { rng_initialize(0); }
};
#endif

/************************************************
 * functions
 ************************************************/

/**
 * @brief Register MPI callbacks of thermostat objects
 *
 * @param thermostat        The thermostat object name
 */
#define NEW_THERMOSTAT(thermostat)                                             \
  void mpi_##thermostat##_set_rng_seed(uint32_t seed);                         \
  void thermostat##_set_rng_seed(uint32_t seed);                               \
  void thermostat##_set_rng_counter(uint64_t seed);

NEW_THERMOSTAT(langevin)
NEW_THERMOSTAT(brownian)
#ifdef NPT
NEW_THERMOSTAT(npt_iso)
#endif
NEW_THERMOSTAT(thermalized_bond)
#ifdef DPD
NEW_THERMOSTAT(dpd)
#endif
#ifdef STOKESIAN_DYNAMICS
NEW_THERMOSTAT(stokesian)
#endif

/* Exported thermostat globals */
extern LangevinThermostat langevin;
extern BrownianThermostat brownian;
#ifdef NPT
extern IsotropicNptThermostat npt_iso;
#endif
extern ThermalizedBondThermostat thermalized_bond;
#ifdef DPD
extern DPDThermostat dpd;
#endif
#ifdef STOKESIAN_DYNAMICS
extern StokesianThermostat stokesian;
#endif

/** Initialize constants of the thermostat at the start of integration */
void thermo_init(double time_step);

/** Increment RNG counters */
void philox_counter_increment();

void mpi_set_brownian_gamma(Thermostat::GammaType const &gamma);
void mpi_set_brownian_gamma_rot(Thermostat::GammaType const &gamma);

void mpi_set_langevin_gamma(Thermostat::GammaType const &gamma);
void mpi_set_langevin_gamma_rot(Thermostat::GammaType const &gamma);

void mpi_set_thermo_virtual(bool thermo_virtual);

void mpi_set_temperature(double temperature);
void mpi_set_temperature_local(double temperature);

void mpi_set_thermo_switch(int thermo_switch);
void mpi_set_thermo_switch_local(int thermo_switch);

#ifdef NPT
void mpi_set_nptiso_gammas(double gamma0, double gammav);
void mpi_set_nptiso_gammas_local(double gamma0, double gammav);
#endif

#endif
