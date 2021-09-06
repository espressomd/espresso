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

#include "config.hpp"

#ifdef NPT
#include "velocity_verlet_npt.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "rotation.hpp"
#include "thermostat.hpp"
#include "thermostats/npt_inline.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives.hpp>

#include <cmath>
#include <functional>

void velocity_verlet_npt_propagate_vel_final(const ParticleRange &particles,
                                             double time_step) {
  extern IsotropicNptThermostat npt_iso;
  nptiso.p_vel = {};

  for (auto &p : particles) {
    // Virtual sites are not propagated during integration
    if (p.p.is_virtual)
      continue;
    auto const noise = friction_therm0_nptiso<2>(npt_iso, p.m.v, p.p.identity);
    for (int j = 0; j < 3; j++) {
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
        if (nptiso.geometry & nptiso.nptgeom_dir[j]) {
          nptiso.p_vel[j] += Utils::sqr(p.m.v[j] * time_step) * p.p.mass;
          p.m.v[j] += (p.f.f[j] * time_step / 2.0 + noise[j]) / p.p.mass;
        } else {
          // Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt)
          p.m.v[j] += p.f.f[j] * time_step / 2.0 / p.p.mass;
        }
      }
    }
  }
}

/** Scale and communicate instantaneous NpT pressure */
void velocity_verlet_npt_finalize_p_inst(double time_step) {
  extern IsotropicNptThermostat npt_iso;
  /* finalize derivation of p_inst */
  nptiso.p_inst = 0.0;
  for (int i = 0; i < 3; i++) {
    if (nptiso.geometry & nptiso.nptgeom_dir[i]) {
      nptiso.p_vel[i] /= Utils::sqr(time_step);
      nptiso.p_inst += nptiso.p_vir[i] + nptiso.p_vel[i];
    }
  }

  double p_sum = 0.0;
  boost::mpi::reduce(comm_cart, nptiso.p_inst, p_sum, std::plus<double>(), 0);
  if (this_node == 0) {
    nptiso.p_inst = p_sum / (nptiso.dimension * nptiso.volume);
    nptiso.p_diff += (nptiso.p_inst - nptiso.p_ext) * 0.5 * time_step +
                     friction_thermV_nptiso(npt_iso, nptiso.p_diff);
  }
}

void velocity_verlet_npt_propagate_pos(const ParticleRange &particles,
                                       double time_step) {
  Utils::Vector3d scal{};
  double L_new = 0.0;

  /* finalize derivation of p_inst */
  velocity_verlet_npt_finalize_p_inst(time_step);

  /* adjust \ref NptIsoParameters::nptiso.volume; prepare pos- and
   * vel-rescaling
   */
  if (this_node == 0) {
    nptiso.volume += nptiso.inv_piston * nptiso.p_diff * 0.5 * time_step;
    scal[2] = Utils::sqr(box_geo.length()[nptiso.non_const_dim]) /
              pow(nptiso.volume, 2.0 / nptiso.dimension);
    nptiso.volume += nptiso.inv_piston * nptiso.p_diff * 0.5 * time_step;
    if (nptiso.volume < 0.0) {

      runtimeErrorMsg()
          << "your choice of piston= " << nptiso.piston << ", dt= " << time_step
          << ", p_diff= " << nptiso.p_diff
          << " just caused the volume to become negative, decrease dt";
      nptiso.volume = box_geo.volume();
      scal[2] = 1;
    }

    L_new = pow(nptiso.volume, 1.0 / nptiso.dimension);

    scal[1] = L_new * box_geo.length_inv()[nptiso.non_const_dim];
    scal[0] = 1. / scal[1];
  }
  boost::mpi::broadcast(comm_cart, scal, 0);

  /* propagate positions while rescaling positions and velocities */
  for (auto &p : particles) {
    if (p.p.is_virtual)
      continue;
    for (int j = 0; j < 3; j++) {
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
        if (nptiso.geometry & nptiso.nptgeom_dir[j]) {
          p.r.p[j] = scal[1] * (p.r.p[j] + scal[2] * p.m.v[j] * time_step);
          p.l.p_old[j] *= scal[1];
          p.m.v[j] *= scal[0];
        } else {
          p.r.p[j] += p.m.v[j] * time_step;
        }
      }
    }
  }

  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);

  /* Apply new volume to the box-length, communicate it, and account for
   * necessary adjustments to the cell geometry */
  Utils::Vector3d new_box;

  if (this_node == 0) {
    new_box = box_geo.length();

    for (int i = 0; i < 3; i++) {
      if (nptiso.cubic_box || nptiso.geometry & nptiso.nptgeom_dir[i]) {
        new_box[i] = L_new;
      }
    }
  }

  boost::mpi::broadcast(comm_cart, new_box, 0);

  box_geo.set_length(new_box);
  // fast box length update
  on_boxl_change(true);
}

void velocity_verlet_npt_propagate_vel(const ParticleRange &particles,
                                       double time_step) {
  extern IsotropicNptThermostat npt_iso;
  nptiso.p_vel = {};

  for (auto &p : particles) {
#ifdef ROTATION
    propagate_omega_quat_particle(p, time_step);
#endif

    // Don't propagate translational degrees of freedom of vs
    if (p.p.is_virtual)
      continue;
    for (int j = 0; j < 3; j++) {
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
        auto const noise =
            friction_therm0_nptiso<1>(npt_iso, p.m.v, p.p.identity);
        if (integ_switch == INTEG_METHOD_NPT_ISO &&
            (nptiso.geometry & nptiso.nptgeom_dir[j])) {
          p.m.v[j] += (p.f.f[j] * time_step / 2.0 + noise[j]) / p.p.mass;
          nptiso.p_vel[j] += Utils::sqr(p.m.v[j] * time_step) * p.p.mass;
        } else {
          // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * a(t)
          p.m.v[j] += p.f.f[j] * time_step / 2.0 / p.p.mass;
        }
      }
    }
  }
}

void velocity_verlet_npt_step_1(const ParticleRange &particles,
                                double time_step) {
  velocity_verlet_npt_propagate_vel(particles, time_step);
  velocity_verlet_npt_propagate_pos(particles, time_step);
  increment_sim_time(time_step);
}

void velocity_verlet_npt_step_2(const ParticleRange &particles,
                                double time_step) {
  velocity_verlet_npt_propagate_vel_final(particles, time_step);
#ifdef ROTATION
  convert_torques_propagate_omega(particles, time_step);
#endif
  velocity_verlet_npt_finalize_p_inst(time_step);
}
#endif
