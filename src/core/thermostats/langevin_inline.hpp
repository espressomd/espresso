/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#include "Particle.hpp"
#include "random.hpp"
#include "rotation.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/attributes.hpp>
#include <utils/matrix.hpp>

// The Langevin friction kernels have column-kernel overloads (below) that read
// the velocity / angular-velocity / id (+ gamma, quaternion) columns via
// hoisted *_view() handles indexed by row. The RNG key uses the id column as
// the Philox key1, so the per-particle noise stream is bitwise identical
// regardless of iteration order. The Particle-view overloads are retained for
// non-hot / unit-test callers.

// Force-inline the Langevin friction kernels: friction_thermo_langevin must
// stay within gcc's inline budget for the init_forces_and_thermostat lambda
// (forces.cpp) to avoid an out-of-line call overhead. See
// src/utils/include/utils/attributes.hpp for the macro definition.

/** Langevin thermostat for particle translational velocities.
 *  @param[in]     langevin       Parameters
 *  @param[in]     p              Particle
 *  @param[in]     time_step      Time step
 *  @param[in]     kT             Thermal energy
 */
ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3d
friction_thermo_langevin(LangevinThermostat const &langevin, Particle const &p,
                         double time_step, double kT) {
  using namespace Thermostat;
  // Determine prefactors for the friction and the noise term
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto const gamma = handle_particle_gamma(p.gamma(), langevin.gamma);
  auto const pref_friction = -gamma;
  auto const pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
#else
  auto const pref_friction = langevin.pref_friction;
  auto const pref_noise = langevin.pref_noise;
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN>(
      langevin.rng_counter(), langevin.rng_seed(), p.id());
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  auto const aniso_flag = (pref_friction[0] != pref_friction[1]) ||
                          (pref_friction[1] != pref_friction[2]);
  if (aniso_flag) {
    // Apply (O * diag(g) * O^T) * v as O * (g ⊙ (O^T * v)) using a single
    // rotation matrix shared between the friction and noise terms.
    auto const O = rotation_matrix(Utils::Quaternion<double>(p.quat()));
    auto const Ot = O.transposed();
    auto const v_body = Ot * Utils::Vector3d(p.v());
    auto const noise_body = Ot * noise;
    return O * (hadamard_product(pref_friction, v_body) +
                hadamard_product(pref_noise, noise_body));
  }
#endif
  return hadamard_product(pref_friction, Utils::Vector3d(p.v())) +
         hadamard_product(pref_noise, noise);
}

#ifdef ESPRESSO_ROTATION
/** Langevin thermostat for particle angular velocities.
 *  @param[in]     langevin       Parameters
 *  @param[in]     p              Particle
 *  @param[in]     time_step      Time step
 *  @param[in]     kT             Thermal energy
 */
ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3d
friction_thermo_langevin_rotation(LangevinThermostat const &langevin,
                                  Particle const &p, double time_step,
                                  double kT) {
  using namespace Thermostat;

#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto const gamma =
      handle_particle_gamma(p.gamma_rot(), langevin.gamma_rotation);
  auto const pref_friction = gamma;
  auto const pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
#else
  auto const pref_friction = langevin.gamma_rotation;
  auto const pref_noise = langevin.pref_noise_rotation;
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(
      langevin.rng_counter(), langevin.rng_seed(), p.id());
  return -hadamard_product(pref_friction, Utils::Vector3d(p.omega())) +
         hadamard_product(pref_noise, noise);
}
#endif // ESPRESSO_ROTATION

// -- column-kernel overloads ----------------------------------------------
// Same arithmetic as the Particle-view overloads above; velocity / omega / id
// (+ gamma, quaternion) are read from hoisted *_view() handles by row. The RNG
// key uses id_view(row) as the Philox key1.

/** Langevin thermostat for particle translational velocities (column form). */
template <class VelView, class IdView
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
          ,
          class GammaView
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
          ,
          class QuatView
#endif
          >
ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3d
friction_thermo_langevin(LangevinThermostat const &langevin,
                         VelView const &vel_view, IdView const &id_view,
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                         GammaView const &gamma_view,
#endif
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
                         QuatView const &quat_view,
#endif
                         int const row, double time_step, double kT) {
  using namespace Thermostat;
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  ParticleStore::GammaValue const p_gamma{
      gamma_view(row, 0), gamma_view(row, 1), gamma_view(row, 2)};
#else
  ParticleStore::GammaValue const p_gamma = gamma_view(row);
#endif
  auto const gamma = handle_particle_gamma(p_gamma, langevin.gamma);
  auto const pref_friction = -gamma;
  auto const pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
#else
  auto const pref_friction = langevin.pref_friction;
  auto const pref_noise = langevin.pref_noise;
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN>(
      langevin.rng_counter(), langevin.rng_seed(), id_view(row));
  Utils::Vector3d const v{vel_view(row, 0), vel_view(row, 1), vel_view(row, 2)};
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  auto const aniso_flag = (pref_friction[0] != pref_friction[1]) ||
                          (pref_friction[1] != pref_friction[2]);
  if (aniso_flag) {
    Utils::Quaternion<double> quat;
    for (unsigned int j = 0u; j < 4u; ++j)
      quat[j] = quat_view(row, j);
    auto const O = rotation_matrix(quat);
    auto const Ot = O.transposed();
    auto const v_body = Ot * v;
    auto const noise_body = Ot * noise;
    return O * (hadamard_product(pref_friction, v_body) +
                hadamard_product(pref_noise, noise_body));
  }
#endif
  return hadamard_product(pref_friction, v) +
         hadamard_product(pref_noise, noise);
}

#ifdef ESPRESSO_ROTATION
/** Langevin thermostat for particle angular velocities (column form). */
template <class OmegaView, class IdView
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
          ,
          class GammaRotView
#endif
          >
ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3d
friction_thermo_langevin_rotation(LangevinThermostat const &langevin,
                                  OmegaView const &omega_view,
                                  IdView const &id_view,
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
                                  GammaRotView const &gamma_rot_view,
#endif
                                  int const row, double time_step, double kT) {
  using namespace Thermostat;
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  ParticleStore::GammaValue const p_gamma_rot{
      gamma_rot_view(row, 0), gamma_rot_view(row, 1), gamma_rot_view(row, 2)};
#else
  ParticleStore::GammaValue const p_gamma_rot = gamma_rot_view(row);
#endif
  auto const gamma =
      handle_particle_gamma(p_gamma_rot, langevin.gamma_rotation);
  auto const pref_friction = gamma;
  auto const pref_noise = LangevinThermostat::sigma(kT, time_step, gamma);
#else
  auto const pref_friction = langevin.gamma_rotation;
  auto const pref_noise = langevin.pref_noise_rotation;
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

  auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(
      langevin.rng_counter(), langevin.rng_seed(), id_view(row));
  Utils::Vector3d const omega{omega_view(row, 0), omega_view(row, 1),
                              omega_view(row, 2)};
  return -hadamard_product(pref_friction, omega) +
         hadamard_product(pref_noise, noise);
}
#endif // ESPRESSO_ROTATION
