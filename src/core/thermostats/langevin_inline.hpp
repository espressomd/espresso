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
  Thermostat::GammaType pref_friction = langevin.pref_friction;
  Thermostat::GammaType pref_noise = langevin.pref_noise;
#ifdef THERMOSTAT_PER_PARTICLE
  // override default if particle-specific gamma
  if (p.gamma() >= Thermostat::GammaType{}) {
    auto const gamma =
        p.gamma() >= Thermostat::GammaType{} ? p.gamma() : langevin.gamma;
    pref_friction = -gamma;
    pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
  }
#endif // THERMOSTAT_PER_PARTICLE

  // Get effective velocity in the thermostatting
#ifdef ENGINE
  auto const &velocity = (p.swimming().v_swim != 0)
                             ? p.v() - p.swimming().v_swim * p.calc_director()
                             : p.v();
#else
  auto const &velocity = p.v();
#endif // ENGINE
#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check
  auto const aniso_flag = (pref_friction[0] != pref_friction[1]) ||
                          (pref_friction[1] != pref_friction[2]);

  // In case of anisotropic particle: body-fixed reference frame. Otherwise:
  // lab-fixed reference frame.
  const Utils::Matrix<double, 3, 3> fric_mat =
      boost::qvm::diag_mat(pref_friction);
  const Utils::Matrix<double, 3, 3> noise_mat =
      boost::qvm::diag_mat(pref_noise);

  auto const friction_op =
      aniso_flag ? convert_body_to_space(p, fric_mat) : fric_mat;
  auto const noise_op =
      aniso_flag ? convert_body_to_space(p, noise_mat) : noise_mat;
#else
  auto const &friction_op = pref_friction;
  auto const &noise_op = pref_noise;
#endif // PARTICLE_ANISOTROPY

  return friction_op * velocity +
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

  auto pref_friction = -langevin.gamma_rotation;
  auto pref_noise = langevin.pref_noise_rotation;

#ifdef THERMOSTAT_PER_PARTICLE
  // override default if particle-specific gamma
  if (p.gamma_rot() >= Thermostat::GammaType{}) {
    auto const gamma = p.gamma_rot() >= Thermostat::GammaType{}
                           ? p.gamma_rot()
                           : langevin.gamma_rotation;
    pref_friction = -gamma;
    pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
  }
#endif // THERMOSTAT_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(
      langevin.rng_counter(), langevin.rng_seed(), p.id());
  return hadamard_product(pref_friction, p.omega()) +
         hadamard_product(pref_noise, noise);
}

#endif // ROTATION
#endif
