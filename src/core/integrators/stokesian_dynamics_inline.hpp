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

#ifdef STOKESIAN_DYNAMICS

#include "rotation.hpp"
#include "stokesian_dynamics/sd_interface.hpp"
#include "thermostat.hpp"

inline void stokesian_dynamics_step_1(ParticleRangeStokesian const &particles,
                                      StokesianThermostat const &stokesian,
                                      double time_step, double kT) {
  propagate_vel_pos_sd(particles, stokesian, time_step, kT);

  for (auto &p : particles) {
    // translate
    p.pos() += p.v() * time_step;
    // rotate
    auto const norm = p.omega().norm();
    if (norm != 0.) {
      auto const omega_unit = (1. / norm) * p.omega();
      local_rotate_particle(p, omega_unit, norm * time_step);
    }
  }
}

#endif // STOKESIAN_DYNAMICS
