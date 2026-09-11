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

#include "cell_system/CellStructure.hpp"
#include "system/System.hpp"

#include <vector>

bool add_bond(System::System &system, int bond_id,
              std::vector<int> const &particle_ids) {
  // Mutable optional view -- the bond is inserted through it below.
  auto p = system.cell_structure->get_local_particle(particle_ids[0]);
  if (p) {
    // The bond view is stored in the bond list of the primary particle.
    // Thus the bond views' partner list only contains the other particle id.
    BondView bond(bond_id, {particle_ids.data() + 1, particle_ids.size() - 1});
    p->bonds().insert(bond);
    return true;
  }
  return false;
}
