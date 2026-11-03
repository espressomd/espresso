/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

// Non-template wrapper around update_cabana_state<VerletCriterion<>>. The
// AoSoA-commit / Verlet-list-build template giant is instantiated here ONCE
// (in short_range_verlet.cpp) instead of separately in forces.cpp, energy.cpp
// and pressure.cpp. Those TUs call this wrapper, so the giant no longer lands
// in each of their objects.

namespace System {
class System;
}

// Builds the Verlet criterion from the system's cutoffs and updates the
// cabana state. Only the collision-detection cutoff is caller-specific
// (energy and pressure observables pass the inactive cutoff); everything
// else is derived from @p system so the call sites cannot drift. When no
// electrostatics/dipolar/collision cutoff is active, the pure-short-range
// criterion variant is used, dropping those dead branches from the
// per-candidate build loop.
void update_verlet_state(System::System const &system,
                         double collision_detection_cutoff);
