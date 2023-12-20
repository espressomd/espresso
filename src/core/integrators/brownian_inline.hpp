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

#pragma once

#include "config/config.hpp"

#include "rotation.hpp"
#include "thermostat.hpp"
#include "thermostats/brownian_inline.hpp"

inline void brownian_dynamics_propagator(BrownianThermostat const &brownian,
                                         Particle &p, double time_step,
                                         double kT) {
  p.pos() += bd_drag(brownian.gamma, p, time_step);
  p.v() = bd_drag_vel(brownian.gamma, p);
  p.pos() += bd_random_walk(brownian, p, time_step, kT);
  p.v() += bd_random_walk_vel(brownian, p);
}

#ifdef ROTATION
inline void brownian_dynamics_rotator(BrownianThermostat const &brownian,
                                      Particle &p, double time_step,
                                      double kT) {
  if (!p.can_rotate())
    return;
  convert_torque_to_body_frame_apply_fix(p);
  p.quat() = bd_drag_rot(brownian.gamma_rotation, p, time_step);
  p.omega() = bd_drag_vel_rot(brownian.gamma_rotation, p);
  p.quat() = bd_random_walk_rot(brownian, p, time_step, kT);
  p.omega() += bd_random_walk_vel_rot(brownian, p);
}
#endif // ROTATION
