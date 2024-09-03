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

class GlueToSurface {
public:
  double distance;
  // Square of distance at which particle are bound
  double distance_sq;
  /// bond type used between centers of colliding particles
  int bond_centers;
  /// bond type used between virtual sites
  int bond_vs;
  /// particle type for virtual sites created on collision
  int part_type_vs;
  /// For mode "glue to surface": The distance from the particle which is to be
  /// glued to the new virtual site
  double dist_glued_part_to_vs;
  /// For mode "glue to surface": The particle type being glued
  int part_type_to_be_glued;
  /// For mode "glue to surface": The particle type to which the virtual site is
  /// attached
  int part_type_to_attach_vs_to;
  /// Particle type to which the newly glued particle is converted
  int part_type_after_glueing;

  GlueToSurface(double distance, int bond_centers, int bond_vs,
                int part_type_vs, double dist_glued_part_to_vs,
                int part_type_to_be_glued, int part_type_to_attach_vs_to,
                int part_type_after_glueing)
      : distance{distance}, distance_sq{distance * distance},
        bond_centers{bond_centers}, bond_vs{bond_vs},
        part_type_vs{part_type_vs},
        dist_glued_part_to_vs{dist_glued_part_to_vs},
        part_type_to_be_glued{part_type_to_be_glued},
        part_type_to_attach_vs_to{part_type_to_attach_vs_to},
        part_type_after_glueing{part_type_after_glueing} {}

  void initialize(System::System &system);

  auto cutoff() const { return distance; }

  /** @brief Check additional criteria for the glue_to_surface collision mode */
  bool glue_to_surface_criterion(Particle const &p1, Particle const &p2) const {
    return (((p1.type() == part_type_to_be_glued) and
             (p2.type() == part_type_to_attach_vs_to)) or
            ((p2.type() == part_type_to_be_glued) and
             (p1.type() == part_type_to_attach_vs_to)));
  }

  void handle_collisions(System::System &system,
                         std::vector<CollisionPair> &local_collision_queue);

  bool detect_collision(Particle const &p1, Particle const &p2,
                        double const dist_sq) const {
    if (dist_sq > distance_sq)
      return false;

    if (!glue_to_surface_criterion(p1, p2))
      return false;

    return detect_collision_common(p1, p2, bond_centers);
  }
};

} // namespace CollisionDetection
#endif // VIRTUAL_SITES_RELATIVE
#endif // COLLISION_DETECTION
