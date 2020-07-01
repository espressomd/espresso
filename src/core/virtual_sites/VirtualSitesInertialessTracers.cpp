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
#include "grid_based_algorithms/lb_interface.hpp"
#include "virtual_sites/lb_inertialess_tracers.hpp"
#include <algorithm>

void VirtualSitesInertialessTracers::after_force_calc() {
  // Now the forces are computed and need to go into the LB fluid
  if (lattice_switch == ActiveLB::CPU) {
    IBM_ForcesIntoFluid_CPU();
    return;
  }
#ifdef CUDA
  if (lattice_switch == ActiveLB::GPU) {
    IBM_ForcesIntoFluid_GPU(cell_structure.local_particles());
    return;
  }
#endif
  if (std::any_of(cell_structure.local_particles().begin(),
                  cell_structure.local_particles().end(),
                  [](Particle &p) { return p.p.is_virtual; })) {
    runtimeErrorMsg() << "Inertialess Tracers: No LB method was active but "
                         "virtual sites present.";
    return;
  }
}

void VirtualSitesInertialessTracers::after_lb_propagation() {
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  IBM_UpdateParticlePositions(cell_structure.local_particles());
#endif // VS inertialess tracers
}
#endif
