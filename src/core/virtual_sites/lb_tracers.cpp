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
#include "config/config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "BoxGeometry.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "lb/particle_coupling.hpp"
#include "system/System.hpp"

struct DeferredActiveLBChecks {
  DeferredActiveLBChecks() {
    auto const &lb = System::get_system().lb;
    m_value = lb.is_solver_set();
  }
  auto operator()() const { return m_value; }
  bool m_value;
};

static bool lb_active_check(DeferredActiveLBChecks const &check) {
  if (not check()) {
    runtimeErrorMsg() << "LB needs to be active for inertialess tracers.";
    return false;
  }
  return true;
}

void lb_tracers_add_particle_force_to_fluid(CellStructure &cell_structure,
                                            double time_step) {
  DeferredActiveLBChecks check_lb_solver_set{};
  auto &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto const agrid = (check_lb_solver_set()) ? system.lb.get_agrid() : 0.;
  auto const to_lb_units = (check_lb_solver_set()) ? 1. / agrid : 0.;

  // Distribute summed-up forces from physical particles to ghosts
  init_forces_ghosts(cell_structure.ghost_particles());
  cell_structure.update_ghosts_and_resort_particle(Cells::DATA_PART_FORCE);

  // Keep track of ghost particles (ids) that have already been coupled
  LB::CouplingBookkeeping bookkeeping{};
  // Apply particle forces to the LB fluid at particle positions.
  // For physical particles, also set particle velocity = fluid velocity.
  for (auto const &particle_range :
       {cell_structure.local_particles(), cell_structure.ghost_particles()}) {
    for (auto const &p : particle_range) {
      if (!LB::is_tracer(p))
        continue;
      if (!lb_active_check(check_lb_solver_set)) {
        return;
      }
      if (bookkeeping.should_be_coupled(p)) {
        for (auto pos : positions_in_halo(p.pos(), box_geo, agrid)) {
          add_md_force(system.lb, pos * to_lb_units, p.force(), time_step);
        }
      }
    }
  }

  // Clear ghost forces to avoid double counting later
  init_forces_ghosts(cell_structure.ghost_particles());
}

void lb_tracers_propagate(CellStructure &cell_structure, double time_step) {
  DeferredActiveLBChecks check_lb_solver_set{};
  auto &system = System::get_system();
  auto const &lb = system.lb;
  auto const verlet_skin = cell_structure.get_verlet_skin();
  auto const verlet_skin_sq = verlet_skin * verlet_skin;

  // Advect particles
  for (auto &p : cell_structure.local_particles()) {
    if (!LB::is_tracer(p))
      continue;
    if (!lb_active_check(check_lb_solver_set)) {
      return;
    }
    p.v() = lb.get_coupling_interpolated_velocity(p.pos());
    for (unsigned int i = 0u; i < 3u; i++) {
      if (!p.is_fixed_along(i)) {
        p.pos()[i] += p.v()[i] * time_step;
      }
    }
    // Verlet list update check
    if ((p.pos() - p.pos_at_last_verlet_update()).norm2() > verlet_skin_sq) {
      cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
    }
  }
}
#endif // VIRTUAL_SITES_INERTIALESS_TRACERS
