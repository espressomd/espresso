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

#include "CollisionPair.hpp"
#include "common.hpp"

#include "Particle.hpp"
#include "system/System.hpp"

#include <vector>

namespace CollisionDetection {

class BindCenters {
public:
  /// Distance at which particle are bound
  double distance;
  /// Square of distance at which particle are bound
  double distance_sq;
  /// bond type used between centers of colliding particles
  int bond_centers;

  BindCenters(double distance, int bond_centers)
      : distance{distance}, distance_sq{distance * distance},
        bond_centers{bond_centers} {}

  void initialize(System::System &system);

  auto cutoff() const { return distance; }

  void handle_collisions(System::System &system,
                         std::vector<CollisionPair> &local_collision_queue);

  bool detect_collision(Particle const &p1, Particle const &p2,
                        double const dist_sq) const {
    if (dist_sq > distance_sq)
      return false;

    return detect_collision_common(p1, p2, bond_centers);
  }
};

} // namespace CollisionDetection
#endif // COLLISION_DETECTION
