/*
 * Copyright (C) 2011-2024 The ESPResSo project
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

#include <config/config.hpp>

#ifdef COLLISION_DETECTION

#include "BindCenters.hpp"
#include "CollisionPair.hpp"
#include "utils.hpp"

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "system/System.hpp"

#include <utils/math/sqr.hpp>

#include <cassert>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CollisionDetection {

void BindCenters::initialize(System::System &system) {
  // Validate distance
  if (distance <= 0.) {
    throw std::domain_error("Parameter 'distance' must be > 0");
  }
  // Cache square of cutoff
  distance_sq = Utils::sqr(distance);
  // Check if bond exists
  assert(system.bonded_ias->contains(bond_centers));
  // If the bond type to bind particle centers is not a pair bond...
  // Check that the bonds have the right number of partners
  if (number_of_partners(*system.bonded_ias->at(bond_centers)) != 1) {
    throw std::runtime_error("The bond type to be used for binding particle "
                             "centers needs to be a pair bond");
  }
}

void BindCenters::handle_collisions(
    System::System &system, std::vector<CollisionPair> &local_collision_queue) {
  add_bind_centers(local_collision_queue, *system.cell_structure, bond_centers);
}

} // namespace CollisionDetection

#endif // COLLISION_DETECTION
