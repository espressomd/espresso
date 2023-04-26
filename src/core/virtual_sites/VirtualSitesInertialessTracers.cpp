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

#include "VirtualSitesInertialessTracers.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "virtual_sites/lb_inertialess_tracers.hpp"

#include <algorithm>

static void check_no_vs_exist(char const *const message) {
  if (std::any_of(cell_structure.local_particles().begin(),
                  cell_structure.local_particles().end(),
                  [](Particle const &p) { return p.is_virtual(); })) {
    runtimeErrorMsg() << "Inertialess Tracers: " << message;
  }
}

void VirtualSitesInertialessTracers::after_force_calc() {
  // Now the forces are computed and need to go into the LB fluid
  if (lattice_switch == ActiveLB::CPU) {
    IBM_ForcesIntoFluid_CPU();
    return;
  }
#ifdef CUDA
  if (lattice_switch == ActiveLB::GPU) {
    IBM_ForcesIntoFluid_GPU(cell_structure.local_particles(), this_node);
    if (comm_cart.size() != 1 and this_node != 0) {
      check_no_vs_exist("The LB GPU method cannot integrate virtual sites when "
                        "more than 1 MPI ranks are used. The particles on MPI "
                        "rank >= 2 are now in an undeterminate state.");
    }
    return;
  }
#endif
  check_no_vs_exist("No LB method was active but virtual sites present.");
}

void VirtualSitesInertialessTracers::after_lb_propagation(double time_step) {
  IBM_UpdateParticlePositions(cell_structure.local_particles(), time_step,
                              this_node);
}
#endif // VIRTUAL_SITES_INERTIALESS_TRACERS
