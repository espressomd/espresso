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

#pragma once

#include <config/config.hpp>

#ifdef COLLISION_DETECTION

#include "BondList.hpp"
#include "Particle.hpp"

namespace CollisionDetection {

inline auto detect_collision_common(Particle const &p1, Particle const &p2,
                                    int const bond_centers) {
  // Ignore virtual particles
  if (p1.is_virtual() or p2.is_virtual())
    return false;

  // Check, if there's already a bond between the particles
  if (pair_bond_exists_on(p1.bonds(), p2.id(), bond_centers) or
      pair_bond_exists_on(p2.bonds(), p1.id(), bond_centers))
    return false;

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  // do not create bond between ghost particles
  return not(p1.is_ghost() and p2.is_ghost());
}

} // namespace CollisionDetection
#endif // COLLISION_DETECTION
