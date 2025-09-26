/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include <config/config.hpp>

#include "Particle.hpp"
#include "rotation.hpp"

/** Propagate the velocities and positions. Integration steps before force
 *  calculation of the Symplectic Euler integrator: <br> \f[ v(t+\Delta t) =
 *  v(t) + \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
 *  v(t+\Delta t) \f]
 */
inline void symplectic_euler_propagator_1(Particle &p, double time_step) {
  for (unsigned int j = 0; j < 3; j++) {
    if (!p.is_fixed_along(j)) {
      /* Propagate velocities: v(t+dt) = v(t) + dt * a(t) */
      p.v()[j] += time_step * p.force()[j] / p.mass();

      /* Propagate positions: p(t + dt) = p(t) + dt * v(t+dt) */
      p.pos()[j] += time_step * p.v()[j];
    }
  }
}

/** Final integration step of the Symplectic Euler integrator
 *  For symplectic Euler, there is no second step as all updates
 *  are done in step 1.
 */
inline void symplectic_euler_propagator_2(Particle &, double) {
  // No second step needed for symplectic Euler
  // All propagation is done in step 1
}

#ifdef ESPRESSO_ROTATION
inline void symplectic_euler_rotator_1(Particle &p, double time_step) {
  if (p.can_rotate()) {
    // For rotation, we also use symplectic Euler scheme
    // Update angular velocity first, then orientation
    convert_torque_propagate_omega(p, time_step);
    propagate_omega_quat_particle(p, time_step);
  }
}

inline void symplectic_euler_rotator_2(Particle &, double) {
  // No second step needed for symplectic Euler rotation
}
#endif // ESPRESSO_ROTATION
