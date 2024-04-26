/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include "ParticleRange.hpp"
#include "cell_system/CellStructure.hpp"
#include "rotation.hpp"

/** Propagate the velocities and positions. Integration steps before force
 *  calculation of the Velocity Verlet integrator: <br> \f[ v(t+0.5 \Delta t) =
 *  v(t) + 0.5 \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
 *  v(t+0.5 \Delta t) \f]
 */
inline void velocity_verlet_propagator_1(Particle &p, double time_step) {
  for (unsigned int j = 0; j < 3; j++) {
    if (!p.is_fixed_along(j)) {
      /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5 * dt * a(t) */
      p.v()[j] += 0.5 * time_step * p.force()[j] / p.mass();

      /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
       * v(t+0.5*dt) */
      p.pos()[j] += time_step * p.v()[j];
    }
  }
}

/** Final integration step of the Velocity Verlet integrator
 *  \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f]
 */
inline void velocity_verlet_propagator_2(Particle &p, double time_step) {
  for (unsigned int j = 0; j < 3; j++) {
    if (!p.is_fixed_along(j)) {
      /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
      p.v()[j] += 0.5 * time_step * p.force()[j] / p.mass();
    }
  }
}

#ifdef ROTATION
inline void velocity_verlet_rotator_1(Particle &p, double time_step) {
  if (p.can_rotate())
    propagate_omega_quat_particle(p, time_step);
}

inline void velocity_verlet_rotator_2(Particle &p, double time_step) {
  if (p.can_rotate())
    convert_torque_propagate_omega(p, time_step);
}
#endif // ROTATION
