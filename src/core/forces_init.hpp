/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

// Combined force initialization and Langevin noise application. Defined in its
// own translation unit (forces_init.cpp) so that its [[gnu::flatten]] Langevin
// call tree does not compete with the hot pair kernel and Verlet-list build in
// forces.cpp for gcc's per-TU inline-growth budget.

namespace System {
class System;
}

void init_forces_and_thermostat(System::System const &system);
