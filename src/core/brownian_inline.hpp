/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file brownian_inline.hpp */

#ifndef BROWNIAN_INLINE_HPP
#define BROWNIAN_INLINE_HPP

#include "thermostat.hpp"

#ifdef BROWNIAN_DYNAMICS
/** Propagate position: viscous drag driven by conservative forces.*/
/*********************************************************/
/** \name bd_drag */
/*********************************************************/
/**(Eq. (14.39) T. Schlick, https://doi.org/10.1007/978-1-4419-6351-2 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
inline void bd_drag(Particle &p, double dt) {
  // The friction tensor Z from the Eq. (14.31) of Schlick2010:
  Thermostat::GammaType local_gamma;

  if (p.p.gamma >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma;
  } else {
    local_gamma = langevin_gamma;
  }

  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // Second (deterministic) term of the Eq. (14.39) of Schlick2010.
      // Only a conservative part of the force is used here
#ifndef PARTICLE_ANISOTROPY
      p.r.p[j] += p.f.f[j] * dt / (local_gamma);
#else
      p.r.p[j] += p.f.f[j] * dt / (local_gamma[j]);
#endif // PARTICLE_ANISOTROPY
    }
  }
}

/** Set the terminal velocity driven by the conservative forces drag.*/
/*********************************************************/
/** \name bd_drag_vel */
/*********************************************************/
/**(Eq. (14.34) T. Schlick, https://doi.org/10.1007/978-1-4419-6351-2 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
inline void bd_drag_vel(Particle &p, double dt) {
  // The friction tensor Z from the eq. (14.31) of Schlick2010:
  Thermostat::GammaType local_gamma;

  if (p.p.gamma >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma;
  } else {
    local_gamma = langevin_gamma;
  }

  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (p.p.ext_flag & COORD_FIXED(j)) {
      p.m.v[j] = 0.0;
    } else
#endif
    {
      // First (deterministic) term of the eq. (14.34) of Schlick2010 taking
      // into account eq. (14.35). Only conservative part of the force is used
      // here NOTE: velocity is assigned here and propagated by thermal part
      // further on top of it
#ifndef PARTICLE_ANISOTROPY
      p.m.v[j] = p.f.f[j] / (local_gamma);
#else
      p.m.v[j] = p.f.f[j] / (local_gamma[j]);
#endif // PARTICLE_ANISOTROPY
    }
  }
}

/** Determine the velocities: random walk part.*/
/*********************************************************/
/** \name bd_random_walk_vel */
/*********************************************************/
/**(Eq. (10.2.16) N. Pottier, https://doi.org/10.1007/s10955-010-0114-6 (2010))
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
inline void bd_random_walk_vel(Particle &p, double dt) {
  // skip the translation thermalizing for virtual sites unless enabled
  extern bool thermo_virtual;
  if (p.p.is_virtual && !thermo_virtual)
    return;
  // Just a square root of kT, see eq. (10.2.17) and comments in 2 paragraphs
  // afterwards, Pottier2010
  extern double brown_sigma_vel;
  // first, set defaults
  double brown_sigma_vel_temp = brown_sigma_vel;

  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  auto const constexpr langevin_temp_coeff = 1.0;
  // Is a particle-specific temperature specified?
  if (p.p.T >= 0.) {
    brown_sigma_vel_temp = sqrt(langevin_temp_coeff * p.p.T);
  } else {
    brown_sigma_vel_temp = brown_sigma_vel;
  }
#endif /* LANGEVIN_PER_PARTICLE */

  Utils::Vector3d noise = v_noise_g(p.p.identity, RNGSalt::BROWNIAN);
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // Random (heat) velocity is added here. It is already initialized in the
      // terminal drag part. See eq. (10.2.16) taking into account eq. (10.2.18)
      // and (10.2.29), Pottier2010. Note, that the Pottier2010 units system
      // (see Eq. (10.1.1) there) has been adapted towards the ESPResSo and the
      // referenced above Schlick2010 one, which is defined by the eq. (14.31)
      // of Schlick2010. A difference is the mass factor to the friction tensor.
      // The noise is Gaussian according to the convention at p. 237 (last
      // paragraph), Pottier2010.
      p.m.v[j] += brown_sigma_vel_temp * noise[j] / sqrt(p.p.mass);
    }
  }
}
#endif // BROWNIAN_DYNAMICS

#endif // BROWNIAN_INLINE_HPP
