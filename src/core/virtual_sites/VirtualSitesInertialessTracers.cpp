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

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
#include "VirtualSitesInertialessTracers.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <algorithm>

void VirtualSitesInertialessTracers::after_force_calc(double time_step) {
  // Distribute summed-up forces from physical particles to ghosts
  init_forces_ghosts(cell_structure.ghost_particles());
  cells_update_ghosts(Cells::DATA_PART_FORCE);

  // Apply particle forces to the LB fluid at particle positions
  // For physical particles, also set particle velocity = fluid velocity
  for (auto &p : cell_structure.local_particles()) {
    if (!p.p.is_virtual)
      continue;
    if (lattice_switch == ActiveLB::NONE) {
      runtimeErrorMsg() << "LB needs to be active for inertialess tracers.";
      return;
    };
    if (in_local_halo(p.r.p)) {
      add_md_force(p.r.p / lb_lbfluid_get_agrid(), -p.f.f, time_step);
    }
  };
  for (auto const &p : cell_structure.ghost_particles()) {
    if (!p.p.is_virtual)
      continue;
    if (lattice_switch == ActiveLB::NONE) {
      runtimeErrorMsg() << "LB needs to be active for inertialess tracers.";
      return;
    };
    if (in_local_halo(p.r.p)) {
      add_md_force(p.r.p / lb_lbfluid_get_agrid(), -p.f.f, time_step);
    }
  }

  // Clear ghost forces to avoid double counting later
  init_forces_ghosts(cell_structure.ghost_particles());
}

void VirtualSitesInertialessTracers::after_lb_propagation(double time_step) {
  auto advect_particle = [time_step](auto &p) {
    if (!p.p.is_virtual)
      return;
    if (lattice_switch == ActiveLB::NONE) {
      runtimeErrorMsg() << "LB needs to be active for inertialess tracers.";
      return;
    };
    p.m.v = lb_lbinterpolation_get_interpolated_velocity(p.r.p) *
            lb_lbfluid_get_lattice_speed();
    for (int i = 0; i < 3; i++) {
      if (!(p.p.ext_flag & COORD_FIXED(i))) {
        p.r.p[i] += p.m.v[i] * time_step;
      }
    };
    // verlet list update check
    if ((p.r.p - p.l.p_old).norm2() > skin * skin) {
      cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
    }
  };

  for (auto &p : cell_structure.local_particles())
    advect_particle(p);
}
#endif
