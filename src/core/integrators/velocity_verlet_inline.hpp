#ifndef INTEGRATORS_VELOCITY_VERLET_HPP
#define INTEGRATORS_VELOCITY_VERLET_HPP

#include "config.hpp"
#include "particle_data.hpp"
#include "ParticleRange.hpp"

inline
void velocity_verlet_propagate_vel_pos(const ParticleRange &particles) {
  assert(integ_switch == INTEG_METHOD_NVT);

  auto const skin2 = Utils::sqr(0.5 * skin);
  for (auto &p : particles) {
#ifdef ROTATION
    propagate_omega_quat_particle(&p);
#endif

// Don't propagate translational degrees of freedom of vs
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
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
      set_resort_particles(Cells::RESORT_LOCAL);
  }
}



inline
void velocity_verlet_propagate_vel_final(const ParticleRange &particles) {
  assert(integ_switch == INTEG_METHOD_NVT);

  for (auto &p : particles) {
#ifdef VIRTUAL_SITES
    // Virtual sites are not propagated during integration
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
#endif
        /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
        p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
#ifdef EXTERNAL_FORCES
      }
#endif
    }
  }
}

void velocity_verlet_step_1(const ParticleRange& particles) {
  velocity_verlet_propagate_vel_pos(particles);
  sim_time += time_step;
}

void velocity_verlet_step_2(const ParticleRange& particles) {
  velocity_verlet_propagate_vel_final(particles);
}

#endif
