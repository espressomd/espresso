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

#ifndef THERMOSTATS_LANGEVIN_INLINE_HPP
#define THERMOSTATS_LANGEVIN_INLINE_HPP

#include "config/config.hpp"

#include "Particle.hpp"
#include "random.hpp"
#include "rotation.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

/** Langevin thermostat for particle translational velocities.
 *  Collects the particle velocity (different for ENGINE, PARTICLE_ANISOTROPY).
 *  Collects the langevin parameters kT, gamma (different for
 *  THERMOSTAT_PER_PARTICLE). Applies the noise and friction term.
 *  @param[in]     langevin       Parameters
 *  @param[in]     p              %Particle
 *  @param[in]     time_step      Time step
 *  @param[in]     kT             Temperature
 */
inline Utils::Vector3d
friction_thermo_langevin(LangevinThermostat const &langevin, Particle const &p,
                         double time_step, double kT) {
  // Early exit for virtual particles without thermostat
  if (p.is_virtual() and !thermo_virtual) {
    return {};
  }

  // Determine prefactors for the friction and the noise term
#ifdef THERMOSTAT_PER_PARTICLE
  auto const gamma = handle_particle_gamma(p.gamma(), langevin.gamma);
  auto const pref_friction = -gamma;
  auto const pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
#else
  auto const pref_friction = langevin.pref_friction;
  auto const pref_noise = langevin.pref_noise;
#endif // THERMOSTAT_PER_PARTICLE

  auto const friction_op = handle_particle_anisotropy(p, pref_friction);
  auto const noise_op = handle_particle_anisotropy(p, pref_noise);

  return friction_op * p.v() +
         noise_op * Random::noise_uniform<RNGSalt::LANGEVIN>(
                        langevin.rng_counter(), langevin.rng_seed(), p.id());
}

#ifdef ROTATION
/** Langevin thermostat for particle angular velocities.
 *  Collects the particle velocity (different for PARTICLE_ANISOTROPY).
 *  Collects the langevin parameters kT, gamma_rot (different for
 *  THERMOSTAT_PER_PARTICLE). Applies the noise and friction term.
 *  @param[in]     langevin       Parameters
 *  @param[in]     p              %Particle
 *  @param[in]     time_step      Time step
 *  @param[in]     kT             Temperature
 */
inline Utils::Vector3d
friction_thermo_langevin_rotation(LangevinThermostat const &langevin,
                                  Particle const &p, double time_step,
                                  double kT) {

#ifdef THERMOSTAT_PER_PARTICLE
  auto const gamma =
      handle_particle_gamma(p.gamma_rot(), langevin.gamma_rotation);
  auto const pref_friction = gamma;
  auto const pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
#else
  auto const pref_friction = langevin.gamma_rotation;
  auto const pref_noise = langevin.pref_noise_rotation;
#endif // THERMOSTAT_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(
      langevin.rng_counter(), langevin.rng_seed(), p.id());
  return -hadamard_product(pref_friction, p.omega()) +
         hadamard_product(pref_noise, noise);
}

#endif // ROTATION
#endif
