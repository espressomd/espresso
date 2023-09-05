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
#ifndef STOKESIAN_DYNAMICS_INLINE_HPP
#define STOKESIAN_DYNAMICS_INLINE_HPP

#include "config/config.hpp"

#ifdef STOKESIAN_DYNAMICS
#include "ParticleRange.hpp"
#include "communication.hpp"
#include "integrate.hpp"
#include "rotation.hpp"
#include "stokesian_dynamics/sd_interface.hpp"

inline void stokesian_dynamics_propagate_vel_pos(const ParticleRange &particles,
                                                 double time_step,
                                                 int default_propagation) {

  // Compute new (translational and rotational) velocities
  propagate_vel_pos_sd(particles, comm_cart, time_step, default_propagation);
  int modes = PropagationMode::TRANS_STOKESIAN;
  if (default_propagation & PropagationMode::TRANS_STOKESIAN)
    modes += PropagationMode::TRANS_SYSTEM_DEFAULT;

  for (auto &p : particles) {
    if (!(p.propagation() & modes))
      continue;

    // Translate particle
    p.pos() += p.v() * time_step;

    // Perform rotation
    auto const norm = p.omega().norm();
    if (norm != 0.) {
      auto const omega_unit = (1. / norm) * p.omega();
      local_rotate_particle(p, omega_unit, norm * time_step);
    }
  }
}

inline void stokesian_dynamics_step_1(const ParticleRange &particles,
                                      double time_step,
                                      int default_propagation) {
  stokesian_dynamics_propagate_vel_pos(particles, time_step,
                                       default_propagation);
}

#endif // STOKESIAN_DYNAMICS
#endif
