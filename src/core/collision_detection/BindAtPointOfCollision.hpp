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
#ifdef VIRTUAL_SITES_RELATIVE

#include "CollisionPair.hpp"
#include "common.hpp"

#include "Particle.hpp"
#include "system/System.hpp"

#include <vector>

namespace CollisionDetection {

class BindAtPointOfCollision {
public:
  /// Distance at which particle are bound
  double distance;
  /// Square of distance at which particle are bound
  double distance_sq;
  /// bond type used between centers of colliding particles
  int bond_centers;
  /// bond type used between virtual sites
  int bond_vs;
  /**
   * @brief Barycenter of the virtual site.
   * 0=on same particle as related to,
   * 1=on collision partner,
   * 0.5=in the middle between the two particles
   */
  double vs_placement;
  /// particle type for virtual sites created on collision
  int part_type_vs;

  BindAtPointOfCollision(double distance, int bond_centers, int bond_vs,
                         double vs_placement, int part_type_vs)
      : distance{distance}, distance_sq{distance * distance},
        bond_centers{bond_centers}, bond_vs{bond_vs},
        vs_placement{vs_placement}, part_type_vs{part_type_vs} {}

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
#endif // VIRTUAL_SITES_RELATIVE
#endif // COLLISION_DETECTION
