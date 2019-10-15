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
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"
#include "utils/math/sqr.hpp"

void velocity_verlet_npt_propagate_vel_final(const ParticleRange &particles) {
  nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;

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
        if (nptiso.geometry & nptiso.nptgeom_dir[j]) {
          nptiso.p_vel[j] += Utils::sqr(p.m.v[j] * time_step) * p.p.mass;
          p.m.v[j] += 0.5 * time_step / p.p.mass * p.f.f[j] +
                      friction_therm0_nptiso(p.m.v[j]) / p.p.mass;
        } else
          // Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt)
          p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
#ifdef EXTERNAL_FORCES
      }
#endif
    }
  }
}

/** Scale and communicate instantaneous NpT pressure */
void velocity_verlet_npt_finalize_p_inst() {
  double p_tmp = 0.0;
  int i;
  /* finalize derivation of p_inst */
  nptiso.p_inst = 0.0;
  for (i = 0; i < 3; i++) {
    if (nptiso.geometry & nptiso.nptgeom_dir[i]) {
      nptiso.p_vel[i] /= Utils::sqr(time_step);
      nptiso.p_inst += nptiso.p_vir[i] + nptiso.p_vel[i];
    }
  }

  MPI_Reduce(&nptiso.p_inst, &p_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  if (this_node == 0) {
    nptiso.p_inst = p_tmp / (nptiso.dimension * nptiso.volume);
    nptiso.p_diff = nptiso.p_diff +
                    (nptiso.p_inst - nptiso.p_ext) * 0.5 * time_step +
                    friction_thermV_nptiso(nptiso.p_diff);
  }
}

void velocity_verlet_npt_propagate_pos(const ParticleRange &particles) {
  double scal[3] = {0., 0., 0.}, L_new = 0.0;

  /* finalize derivation of p_inst */
  velocity_verlet_npt_finalize_p_inst();

  /* adjust \ref nptiso_struct::nptiso.volume; prepare pos- and
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
      nptiso.volume =
          box_geo.length()[0] * box_geo.length()[1] * box_geo.length()[2];
      scal[2] = 1;
    }

    L_new = pow(nptiso.volume, 1.0 / nptiso.dimension);

    scal[1] = L_new / box_geo.length()[nptiso.non_const_dim];
    scal[0] = 1 / scal[1];
  }
  MPI_Bcast(scal, 3, MPI_DOUBLE, 0, comm_cart);

  /* propagate positions while rescaling positions and velocities */
  for (auto &p : particles) {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
#endif
        if (nptiso.geometry & nptiso.nptgeom_dir[j]) {
          {
            p.r.p[j] = scal[1] * (p.r.p[j] + scal[2] * p.m.v[j] * time_step);
            p.l.p_old[j] *= scal[1];
            p.m.v[j] *= scal[0];
          }
        } else {
          p.r.p[j] += p.m.v[j] * time_step;
        }
#ifdef EXTERNAL_FORCES
      }
#endif
    }
  }

  set_resort_particles(Cells::RESORT_LOCAL);

  /* Apply new volume to the box-length, communicate it, and account for
   * necessary adjustments to the cell geometry */
  if (this_node == 0) {
    Utils::Vector3d new_box = box_geo.length();

    for (int i = 0; i < 3; i++) {
      if (nptiso.geometry & nptiso.nptgeom_dir[i] || nptiso.cubic_box) {
        new_box[i] = L_new;
      }
    }

    box_geo.set_length(new_box);
  }

  MPI_Bcast(box_geo.m_length.data(), 3, MPI_DOUBLE, 0, comm_cart);

  /* fast box length update */
  grid_changed_box_l(box_geo);
  recalc_maximal_cutoff();
  cells_on_geometry_change(CELL_FLAG_FAST);
}

void velocity_verlet_npt_propagate_vel(const ParticleRange &particles) {
#ifdef NPT
  nptiso.p_vel[0] = nptiso.p_vel[1] = nptiso.p_vel[2] = 0.0;
#endif

  for (auto &p : particles) {
#ifdef ROTATION
    propagate_omega_quat_particle(p);
#endif

#ifdef VIRTUAL_SITES
    // Don't propagate translational degrees of freedom of vs
    if (p.p.is_virtual)
      continue;
#endif
    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
#ifdef NPT
        if (integ_switch == INTEG_METHOD_NPT_ISO &&
            (nptiso.geometry & nptiso.nptgeom_dir[j])) {
          p.m.v[j] += p.f.f[j] * 0.5 * time_step / p.p.mass +
                      friction_therm0_nptiso(p.m.v[j]) / p.p.mass;
          nptiso.p_vel[j] += Utils::sqr(p.m.v[j] * time_step) * p.p.mass;
        } else
#endif
          // Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * a(t)
          p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
      }
    }
  }
}

void velocity_verlet_npt_step_1(const ParticleRange &particles) {
  velocity_verlet_npt_propagate_vel(particles);
  velocity_verlet_npt_propagate_pos(particles);
  sim_time += time_step;
}

void velocity_verlet_npt_step_2(const ParticleRange &particles) {
  velocity_verlet_npt_propagate_vel_final(particles);
#ifdef ROTATION
  convert_torques_propagate_omega(particles);
#endif
  velocity_verlet_npt_finalize_p_inst();
}
#endif
