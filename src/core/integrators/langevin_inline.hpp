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
/** \file */

#ifndef LANGEVIN_INLINE_HPP
#define LANGEVIN_INLINE_HPP

#include "config.hpp"

#include <utils/Vector.hpp>

#include "random.hpp"
#include "thermostat.hpp"

/** Langevin thermostat for particle translational velocities.
 *  Collects the particle velocity (different for ENGINE, PARTICLE_ANISOTROPY).
 *  Collects the langevin parameters kT, gamma (different for
 *  LANGEVIN_PER_PARTICLE). Applies the noise and friction term.
 *  @param[in]     langevin       Parameters
 *  @param[in]     p              %Particle
 */
inline Utils::Vector3d
friction_thermo_langevin(LangevinThermostat const &langevin,
                         Particle const &p) {
  // Early exit for virtual particles without thermostat
  if (p.p.is_virtual && !thermo_virtual) {
    return {};
  }

  // Determine prefactors for the friction and the noise term
  // first, set defaults
  Thermostat::GammaType pref_friction = langevin.pref_friction;
  Thermostat::GammaType pref_noise = langevin.pref_noise;
  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  if (p.p.gamma >= Thermostat::GammaType{} or p.p.T >= 0.) {
    auto const kT = p.p.T >= 0. ? p.p.T : temperature;
    auto const gamma =
        p.p.gamma >= Thermostat::GammaType{} ? p.p.gamma : langevin.gamma;
    pref_friction = -gamma;
    pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
  }
#endif // LANGEVIN_PER_PARTICLE

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
  auto const aniso_flag = (pref_friction[0] != pref_friction[1]) ||
                          (pref_friction[1] != pref_friction[2]);

  // In case of anisotropic particle: body-fixed reference frame. Otherwise:
  // lab-fixed reference frame.
  auto const friction_op =
      aniso_flag ? convert_body_to_space(p, diag_matrix(pref_friction))
                 : diag_matrix(pref_friction);
  auto const noise_op = diag_matrix(pref_noise);
#else
  auto const &friction_op = pref_friction;
  auto const &noise_op = pref_noise;
#endif // PARTICLE_ANISOTROPY

  return friction_op * velocity +
         noise_op * Random::noise_uniform<RNGSalt::LANGEVIN>(langevin.rng_get(),
                                                             p.p.identity);
}

#ifdef ROTATION
/** Langevin thermostat for particle angular velocities.
 *  Collects the particle velocity (different for PARTICLE_ANISOTROPY).
 *  Collects the langevin parameters kT, gamma_rot (different for
 *  LANGEVIN_PER_PARTICLE). Applies the noise and friction term.
 *  @param[in]     langevin       Parameters
 *  @param[in]     p              %Particle
 */
inline Utils::Vector3d
friction_thermo_langevin_rotation(LangevinThermostat const &langevin,
                                  Particle const &p) {

  auto pref_friction = -langevin.gamma_rotation;
  auto pref_noise = langevin.pref_noise_rotation;

  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  if (p.p.gamma_rot >= Thermostat::GammaType{} or p.p.T >= 0.) {
    auto const kT = p.p.T >= 0. ? p.p.T : temperature;
    auto const gamma = p.p.gamma_rot >= Thermostat::GammaType{}
                           ? p.p.gamma_rot
                           : langevin.gamma_rotation;
    pref_friction = -gamma;
    pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
  }
#endif // LANGEVIN_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(
      langevin.rng_get(), p.p.identity);
  return hadamard_product(pref_friction, p.m.omega) +
         hadamard_product(pref_noise, noise);
}

#endif // ROTATION
#endif
