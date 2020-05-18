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

#ifndef BROWNIAN_INLINE_HPP
#define BROWNIAN_INLINE_HPP

#include "config.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include "cells.hpp"
#include "particle_data.hpp"
#include "random.hpp"
#include "thermostat.hpp"

/** Determine position: viscous drag driven by conservative forces.
 *  From eq. (14.39) in @cite Schlick2010.
 *  @param[in]     brownian_gamma Brownian translational gamma
 *  @param[in]     p              %Particle
 *  @param[in]     dt             Time interval
 */
inline Utils::Vector3d bd_drag(Thermostat::GammaType const &brownian_gamma,
                               Particle const &p, double dt) {
  // The friction tensor Z from the Eq. (14.31) of Schlick2010:
  Thermostat::GammaType gamma;

#ifdef BROWNIAN_PER_PARTICLE
  if (p.p.gamma >= Thermostat::GammaType{}) {
    gamma = p.p.gamma;
  } else
#endif
  {
    gamma = brownian_gamma;
  }

#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check.
  auto const aniso_flag = (gamma[0] != gamma[1]) || (gamma[1] != gamma[2]);
  Utils::Vector3d delta_pos_lab;
  if (aniso_flag) {
    auto const force_body = convert_vector_space_to_body(p, p.f.f);
    auto const delta_pos_body = hadamard_division(force_body * dt, gamma);
    delta_pos_lab = convert_vector_body_to_space(p, delta_pos_body);
  }
#endif

  Utils::Vector3d position = {};
  for (int j = 0; j < 3; j++) {
    // Second (deterministic) term of the Eq. (14.39) of Schlick2010.
    // Only a conservative part of the force is used here
#ifdef PARTICLE_ANISOTROPY
    if (aniso_flag) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
        position[j] = delta_pos_lab[j];
      }
    } else {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
        position[j] = p.f.f[j] * dt / gamma[j];
      }
    }
#else
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      position[j] = p.f.f[j] * dt / gamma;
    }
#endif // PARTICLE_ANISOTROPY
  }
  return position;
}

/** Set the terminal velocity driven by the conservative forces drag.
 *  From eq. (14.34) in @cite Schlick2010.
 *  @param[in]     brownian_gamma Brownian translational gamma
 *  @param[in]     p              %Particle
 */
inline Utils::Vector3d bd_drag_vel(Thermostat::GammaType const &brownian_gamma,
                                   Particle const &p) {
  // The friction tensor Z from the eq. (14.31) of Schlick2010:
  Thermostat::GammaType gamma;

#ifdef BROWNIAN_PER_PARTICLE
  if (p.p.gamma >= Thermostat::GammaType{}) {
    gamma = p.p.gamma;
  } else
#endif
  {
    gamma = brownian_gamma;
  }

#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check.
  auto const aniso_flag = (gamma[0] != gamma[1]) || (gamma[1] != gamma[2]);
  Utils::Vector3d vel_lab;
  if (aniso_flag) {
    auto const force_body = convert_vector_space_to_body(p, p.f.f);
    auto const vel_body = hadamard_division(force_body, gamma);
    vel_lab = convert_vector_body_to_space(p, vel_body);
  }
#endif

  Utils::Vector3d velocity = {};
  for (int j = 0; j < 3; j++) {
    // First (deterministic) term of the eq. (14.34) of Schlick2010 taking
    // into account eq. (14.35). Only conservative part of the force is used
    // here.
#ifdef PARTICLE_ANISOTROPY
    if (aniso_flag) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
        velocity[j] = vel_lab[j];
      }
    } else {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
        velocity[j] = p.f.f[j] / gamma[j];
      }
    }
#else // PARTICLE_ANISOTROPY
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      velocity[j] = p.f.f[j] / gamma;
    }
#endif // PARTICLE_ANISOTROPY
  }
  return velocity;
}

/** Determine the positions: random walk part.
 *  From eq. (14.37) in @cite Schlick2010.
 *  @param[in]     brownian       Parameters
 *  @param[in]     p              %Particle
 *  @param[in]     dt             Time interval
 */
inline Utils::Vector3d bd_random_walk(BrownianThermostat const &brownian,
                                      Particle const &p, double dt) {
  // skip the translation thermalizing for virtual sites unless enabled
  extern bool thermo_virtual;
  if (p.p.is_virtual && !thermo_virtual)
    return {};
  // first, set defaults
  Thermostat::GammaType sigma_pos = brownian.sigma_pos;

  // Override defaults if per-particle values for T and gamma are given
#ifdef BROWNIAN_PER_PARTICLE
  if (p.p.gamma >= Thermostat::GammaType{}) {
    // Is a particle-specific temperature also specified?
    if (p.p.T >= 0.) {
      if (p.p.T > 0.0) {
        sigma_pos = BrownianThermostat::sigma(p.p.T, p.p.gamma);
      } else {
        // just an indication of the infinity
        sigma_pos = Thermostat::GammaType{};
      }
    } else
        // default temperature but particle-specific gamma
        if (temperature > 0.0) {
      sigma_pos = BrownianThermostat::sigma(temperature, p.p.gamma);
    } else {
      sigma_pos = Thermostat::GammaType{};
    }
  } // particle-specific gamma
  else {
    // No particle-specific gamma, but is there particle-specific temperature
    if (p.p.T >= 0.) {
      if (p.p.T > 0.0) {
        sigma_pos = BrownianThermostat::sigma(p.p.T, brownian.gamma);
      } else {
        // just an indication of the infinity
        sigma_pos = Thermostat::GammaType{};
      }
    } else {
      // default values for both
      sigma_pos = brownian.sigma_pos;
    }
  }
#endif // BROWNIAN_PER_PARTICLE

  // Eq. (14.37) is factored by the Gaussian noise (12.22) with its squared
  // magnitude defined in the second eq. (14.38), Schlick2010.
  Utils::Vector3d delta_pos_body{};
  auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_WALK>(
      brownian.rng_get(), p.p.identity);
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
#ifndef PARTICLE_ANISOTROPY
      if (sigma_pos > 0.0) {
        delta_pos_body[j] = sigma_pos * sqrt(dt) * noise[j];
      } else {
        delta_pos_body[j] = 0.0;
      }
#else
      if (sigma_pos[j] > 0.0) {
        delta_pos_body[j] = sigma_pos[j] * sqrt(dt) * noise[j];
      } else {
        delta_pos_body[j] = 0.0;
      }
#endif // PARTICLE_ANISOTROPY
    }
  }

#ifdef PARTICLE_ANISOTROPY
  // Particle frictional isotropy check.
  auto const aniso_flag =
      (sigma_pos[0] != sigma_pos[1]) || (sigma_pos[1] != sigma_pos[2]);
  Utils::Vector3d delta_pos_lab;
  if (aniso_flag) {
    delta_pos_lab = convert_vector_body_to_space(p, delta_pos_body);
  }
#endif

  Utils::Vector3d position = {};
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
#ifdef PARTICLE_ANISOTROPY
      position[j] += aniso_flag ? delta_pos_lab[j] : delta_pos_body[j];
#else
      position[j] += delta_pos_body[j];
#endif
    }
  }
  return position;
}

/** Determine the velocities: random walk part.
 *  From eq. (10.2.16) in @cite Pottier2010.
 *  @param[in]     brownian       Parameters
 *  @param[in]     p              %Particle
 */
inline Utils::Vector3d bd_random_walk_vel(BrownianThermostat const &brownian,
                                          Particle const &p) {
  // skip the translation thermalizing for virtual sites unless enabled
  extern bool thermo_virtual;
  if (p.p.is_virtual && !thermo_virtual)
    return {};
  // Just a square root of kT, see eq. (10.2.17) and comments in 2 paragraphs
  // afterwards, Pottier2010
  double sigma_vel;

  // Override defaults if per-particle values for T and gamma are given
#ifdef BROWNIAN_PER_PARTICLE
  // Is a particle-specific temperature specified?
  if (p.p.T >= 0.) {
    sigma_vel = BrownianThermostat::sigma(p.p.T);
  } else {
    sigma_vel = brownian.sigma_vel;
  }
#else
  // defaults
  sigma_vel = brownian.sigma_vel;
#endif // BROWNIAN_PER_PARTICLE

  auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_INC>(
      brownian.rng_get(), p.identity());
  Utils::Vector3d velocity = {};
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // Random (heat) velocity. See eq. (10.2.16) taking into account eq.
      // (10.2.18) and (10.2.29), Pottier2010. Note, that the Pottier2010 units
      // system (see Eq. (10.1.1) there) has been adapted towards the ESPResSo
      // and the referenced above Schlick2010 one, which is defined by the eq.
      // (14.31) of Schlick2010. A difference is the mass factor to the friction
      // tensor. The noise is Gaussian according to the convention at p. 237
      // (last paragraph), Pottier2010.
      velocity[j] += sigma_vel * noise[j] / sqrt(p.p.mass);
    }
  }
  return velocity;
}

#ifdef ROTATION

/** Determine quaternions: viscous drag driven by conservative torques.
 *  An analogy of eq. (14.39) in @cite Schlick2010.
 *  @param[in]     brownian_gamma_rotation Brownian rotational gamma
 *  @param[in]     p              %Particle
 *  @param[in]     dt             Time interval
 */
inline Utils::Vector4d
bd_drag_rot(Thermostat::GammaType const &brownian_gamma_rotation, Particle &p,
            double dt) {
  Thermostat::GammaType gamma;

#ifdef BROWNIAN_PER_PARTICLE
  if (p.p.gamma_rot >= Thermostat::GammaType{}) {
    gamma = p.p.gamma_rot;
  } else
#endif
  {
    gamma = brownian_gamma_rotation;
  }

  Utils::Vector3d dphi = {};
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // only a conservative part of the torque is used here
#ifndef PARTICLE_ANISOTROPY
      dphi[j] = p.f.torque[j] * dt / gamma;
#else
      dphi[j] = p.f.torque[j] * dt / gamma[j];
#endif // ROTATIONAL_INERTIA
    }
  }
  dphi = mask(p.p.rotation, dphi);
  double dphi_m = dphi.norm();
  if (dphi_m != 0.) {
    auto const dphi_u = dphi / dphi_m;
    return local_rotate_particle_body(p, dphi_u, dphi_m);
  }
  return {};
}

/** Set the terminal angular velocity driven by the conservative torques drag.
 *  An analogy of the 1st term of eq. (14.34) in @cite Schlick2010.
 *  @param[in]     brownian_gamma_rotation Brownian rotational gamma
 *  @param[in]     p              %Particle
 */
inline Utils::Vector3d
bd_drag_vel_rot(Thermostat::GammaType const &brownian_gamma_rotation,
                Particle const &p) {
  Thermostat::GammaType gamma;

#ifdef BROWNIAN_PER_PARTICLE
  if (p.p.gamma_rot >= Thermostat::GammaType{}) {
    gamma = p.p.gamma_rot;
  } else
#endif
  {
    gamma = brownian_gamma_rotation;
  }

  Utils::Vector3d omega = {};
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // only conservative part of the force is used here
#ifndef PARTICLE_ANISOTROPY
      omega[j] = p.f.torque[j] / gamma;
#else
      omega[j] = p.f.torque[j] / gamma[j];
#endif // ROTATIONAL_INERTIA
    }
  }
  return mask(p.p.rotation, omega);
}

/** Determine the quaternions: random walk part.
 *  An analogy of eq. (14.37) in @cite Schlick2010.
 *  @param[in]     brownian       Parameters
 *  @param[in]     p              %Particle
 *  @param[in]     dt             Time interval
 */
inline Utils::Vector4d bd_random_walk_rot(BrownianThermostat const &brownian,
                                          Particle const &p, double dt) {
  // first, set defaults
  Thermostat::GammaType sigma_pos = brownian.sigma_pos_rotation;

  // Override defaults if per-particle values for T and gamma are given
#ifdef BROWNIAN_PER_PARTICLE
  if (p.p.gamma_rot >= Thermostat::GammaType{}) {
    // Is a particle-specific temperature also specified?
    if (p.p.T >= 0.) {
      if (p.p.T > 0.0) {
        sigma_pos = BrownianThermostat::sigma(p.p.T, p.p.gamma_rot);
      } else {
        // just an indication of the infinity
        sigma_pos = Utils::Vector3d{};
      }
    } else if (temperature > 0.) {
      // Default temperature but particle-specific gamma
      sigma_pos = BrownianThermostat::sigma(temperature, p.p.gamma_rot);
    } else {
      // just an indication of the infinity
      sigma_pos = Utils::Vector3d{};
    }
  } // particle-specific gamma
  else {
    // No particle-specific gamma, but is there particle-specific temperature
    if (p.p.T >= 0.) {
      if (p.p.T > 0.0) {
        sigma_pos = BrownianThermostat::sigma(p.p.T, brownian.gamma_rotation);
      } else {
        // just an indication of the infinity
        sigma_pos = Utils::Vector3d{};
      }
    } else {
      // Default values for both
      sigma_pos = brownian.sigma_pos_rotation;
    }
  }
#endif // BROWNIAN_PER_PARTICLE

  Utils::Vector3d dphi = {};
  auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_ROT_INC>(
      brownian.rng_get(), p.p.identity);
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
#ifndef PARTICLE_ANISOTROPY
      if (sigma_pos > 0.0) {
        dphi[j] = noise[j] * sigma_pos * sqrt(dt);
      }
#else
      if (sigma_pos[j] > 0.0) {
        dphi[j] = noise[j] * sigma_pos[j] * sqrt(dt);
      }
#endif // ROTATIONAL_INERTIA
    }
  }
  dphi = mask(p.p.rotation, dphi);
  // making the algorithm independent of the order of the rotations
  double dphi_m = dphi.norm();
  if (dphi_m != 0) {
    auto const dphi_u = dphi / dphi_m;
    return local_rotate_particle_body(p, dphi_u, dphi_m);
  }
  return {};
}

/** Determine the angular velocities: random walk part.
 *  An analogy of eq. (10.2.16) in @cite Pottier2010.
 *  @param[in]     brownian       Parameters
 *  @param[in]     p              %Particle
 */
inline Utils::Vector3d
bd_random_walk_vel_rot(BrownianThermostat const &brownian, Particle const &p) {
  double sigma_vel;

  // Override defaults if per-particle values for T and gamma are given
#ifdef BROWNIAN_PER_PARTICLE
  // Is a particle-specific temperature specified?
  if (p.p.T >= 0.) {
    sigma_vel = BrownianThermostat::sigma(p.p.T);
  } else {
    sigma_vel = brownian.sigma_vel_rotation;
  }
#else
  // set defaults
  sigma_vel = brownian.sigma_vel_rotation;
#endif // BROWNIAN_PER_PARTICLE

  Utils::Vector3d domega{};
  auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_ROT_WALK>(
      brownian.rng_get(), p.p.identity);
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      domega[j] = sigma_vel * noise[j] / sqrt(p.p.rinertia[j]);
    }
  }
  return mask(p.p.rotation, domega);
}
#endif // ROTATION

inline void brownian_dynamics_propagator(BrownianThermostat const &brownian,
                                         const ParticleRange &particles) {
  for (auto &p : particles) {
    // Don't propagate translational degrees of freedom of vs
    extern bool thermo_virtual;
    if (!(p.p.is_virtual) or thermo_virtual) {
      p.r.p += bd_drag(brownian.gamma, p, time_step);
      p.m.v = bd_drag_vel(brownian.gamma, p);
      p.r.p += bd_random_walk(brownian, p, time_step);
      p.m.v += bd_random_walk_vel(brownian, p);
      /* Verlet criterion check */
      if ((p.r.p - p.l.p_old).norm2() > Utils::sqr(0.5 * skin))
        cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
    }
#ifdef ROTATION
    if (!p.p.rotation)
      continue;
    convert_torque_to_body_frame_apply_fix(p);
    p.r.quat = bd_drag_rot(brownian.gamma_rotation, p, time_step);
    p.m.omega = bd_drag_vel_rot(brownian.gamma_rotation, p);
    p.r.quat = bd_random_walk_rot(brownian, p, time_step);
    p.m.omega += bd_random_walk_vel_rot(brownian, p);
#endif // ROTATION
  }
  sim_time += time_step;
}

#endif // BROWNIAN_INLINE_HPP
