/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017,2018 The ESPResSo project
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

  if(p.p.gamma >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma;
  } else {
    local_gamma = langevin_gamma;
  }

  double scale_f = 0.5 * time_step * time_step / p.p.mass;
  for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      // Second (deterministic) term of the Eq. (14.39) of Schlick2010.
      // scale_f is required to be aligned with rescaled forces
      // only a conservative part of the force is used here
#ifndef PARTICLE_ANISOTROPY
      p.r.p[j] += p.f.f[j] * dt / (local_gamma * scale_f);
#else
      p.r.p[j] += p.f.f[j] * dt / (local_gamma[j] * scale_f);
#endif // PARTICLE_ANISOTROPY
    }
  }
}
#endif // BROWNIAN_DYNAMICS

#endif // BROWNIAN_INLINE_HPP
