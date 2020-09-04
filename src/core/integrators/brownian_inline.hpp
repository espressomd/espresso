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

#ifndef INTEGRATORS_BROWNIAN_INLINE_HPP
#define INTEGRATORS_BROWNIAN_INLINE_HPP

#include "config.hpp"

#include "CellStructure.hpp"
#include "ParticleRange.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"
#include "thermostats/brownian_inline.hpp"

#include <utils/math/sqr.hpp>

#include <cstdint>

inline void brownian_dynamics_propagator(BrownianThermostat const &brownian,
                                         const ParticleRange &particles,
                                         uint64_t counter) {
  for (auto &p : particles) {
    // Don't propagate translational degrees of freedom of vs
    if (!(p.p.is_virtual) or thermo_virtual) {
      p.r.p += bd_drag(brownian.gamma, p, time_step);
      p.m.v = bd_drag_vel(brownian.gamma, p);
      p.r.p += bd_random_walk(brownian, p, counter, time_step);
      p.m.v += bd_random_walk_vel(brownian, p, counter);
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
    p.r.quat = bd_random_walk_rot(brownian, p, counter, time_step);
    p.m.omega += bd_random_walk_vel_rot(brownian, p, counter);
#endif // ROTATION
  }
  sim_time += time_step;
}

#endif // BROWNIAN_INLINE_HPP
