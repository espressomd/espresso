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
#ifndef INTEGRATORS_VELOCITY_VERLET_HPP
#define INTEGRATORS_VELOCITY_VERLET_HPP

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "config.hpp"
#include "rotation.hpp"

/** Propagate the velocities and positions. Integration steps before force
 *  calculation of the Velocity Verlet integrator: <br> \f[ v(t+0.5 \Delta t) =
 *  v(t) + 0.5 \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
 *  v(t+0.5 \Delta t) \f]
 */
inline void velocity_verlet_propagate_vel_pos(const ParticleRange &particles) {

  auto const skin2 = Utils::sqr(0.5 * skin);
  for (auto &p : particles) {
#ifdef ROTATION
    propagate_omega_quat_particle(p);
#endif

    // Don't propagate translational degrees of freedom of vs
    if (p.p.is_virtual)
      continue;
    for (int j = 0; j < 3; j++) {
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5 * dt * a(t) */
        p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;

        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
         * v(t+0.5*dt) */
        p.r.p[j] += time_step * p.m.v[j];
      }
    }

    /* Verlet criterion check*/
    if (Utils::sqr(p.r.p[0] - p.l.p_old[0]) +
            Utils::sqr(p.r.p[1] - p.l.p_old[1]) +
            Utils::sqr(p.r.p[2] - p.l.p_old[2]) >
        skin2)
      cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

/** Final integration step of the Velocity Verlet integrator
 *  \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f]
 */
inline void
velocity_verlet_propagate_vel_final(const ParticleRange &particles) {

  for (auto &p : particles) {
    // Virtual sites are not propagated during integration
    if (p.p.is_virtual)
      continue;

    for (int j = 0; j < 3; j++) {
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
        /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
        p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
      }
    }
  }
}

inline void velocity_verlet_step_1(const ParticleRange &particles) {
  velocity_verlet_propagate_vel_pos(particles);
  sim_time += time_step;
}

inline void velocity_verlet_step_2(const ParticleRange &particles) {
  velocity_verlet_propagate_vel_final(particles);
#ifdef ROTATION
  convert_torques_propagate_omega(particles);
#endif
}

#endif
