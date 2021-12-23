/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef STOKESIAN_DYNAMICS_INLINE_HPP
#define STOKESIAN_DYNAMICS_INLINE_HPP

#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS
#include "ParticleRange.hpp"
#include "communication.hpp"
#include "integrate.hpp"
#include "rotation.hpp"
#include "stokesian_dynamics/sd_interface.hpp"

inline void stokesian_dynamics_propagate_vel_pos(const ParticleRange &particles,
                                                 double time_step) {

  // Compute new (translational and rotational) velocities
  propagate_vel_pos_sd(particles, comm_cart, time_step);

  for (auto &p : particles) {

    // Translate particle
    p.r.p += p.m.v * time_step;

    // Perform rotation
    auto const norm = p.m.omega.norm();
    if (norm != 0) {
      auto const omega_unit = (1. / norm) * p.m.omega;
      local_rotate_particle(p, omega_unit, norm * time_step);
    }
  }
}

inline void stokesian_dynamics_step_1(const ParticleRange &particles,
                                      double time_step) {
  stokesian_dynamics_propagate_vel_pos(particles, time_step);
  increment_sim_time(time_step);
}

inline void stokesian_dynamics_step_2(const ParticleRange &particles) {}

#endif // STOKESIAN_DYNAMICS
#endif
