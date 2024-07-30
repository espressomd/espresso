/*
 * Copyright (C) 2011-2022 The ESPResSo project
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

#include "config/config.hpp"

#include "cell_system/CellStructure.hpp"

#include "BondList.hpp"
#include "Particle.hpp"

/** @brief Protocols for collision handling. */
enum class CollisionModeType : int {
  /** @brief Deactivate collision detection. */
  OFF = 0,
  /** @brief Create bond between centers of colliding particles. */
  BIND_CENTERS = 1,
  /**
   * @brief Create a bond between the centers of the colliding particles,
   * plus two virtual sites at the point of collision and bind them
   * together. This prevents the particles from sliding against each
   * other. Requires VIRTUAL_SITES_RELATIVE.
   */
  BIND_VS = 2,
  /** @brief Glue a particle to a specific spot on another particle. */
  GLUE_TO_SURF = 3,
};

class Collision_parameters {
public:
  Collision_parameters()
      : mode(CollisionModeType::OFF), distance(0.), distance2(0.),
        bond_centers(-1), bond_vs(-1) {}

  /// collision protocol
  CollisionModeType mode;
  /// distance at which particles are bound
  double distance;
  // Square of distance at which particle are bound
  double distance2;

  /// bond type used between centers of colliding particles
  int bond_centers;
  /// bond type used between virtual sites
  int bond_vs;
  /// particle type for virtual sites created on collision
  int vs_particle_type;

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
  /** Placement of virtual sites for MODE_VS.
   *  0=on same particle as related to,
   *  1=on collision partner,
   *  0.5=in the middle between
   */
  double vs_placement;

  /** @brief Validates parameters and creates particle types if needed. */
  void initialize();
};

/// Parameters for collision detection
extern Collision_parameters collision_params;

#ifdef COLLISION_DETECTION

void prepare_local_collision_queue();

/// Handle the collisions recorded in the queue
void handle_collisions(CellStructure &cell_structure);

/** @brief Add the collision between the given particle ids to the collision
 *  queue
 */
void queue_collision(int part1, int part2);

/** @brief Check additional criteria for the glue_to_surface collision mode */
inline bool glue_to_surface_criterion(Particle const &p1, Particle const &p2) {
  return (((p1.type() == collision_params.part_type_to_be_glued) &&
           (p2.type() == collision_params.part_type_to_attach_vs_to)) ||
          ((p2.type() == collision_params.part_type_to_be_glued) &&
           (p1.type() == collision_params.part_type_to_attach_vs_to)));
}

/** @brief Detect (and queue) a collision between the given particles. */
inline void detect_collision(Particle const &p1, Particle const &p2,
                             double const dist2) {
  if (dist2 > collision_params.distance2)
    return;

  // If we are in the glue to surface mode, check that the particles
  // are of the right type
  if (collision_params.mode == CollisionModeType::GLUE_TO_SURF)
    if (!glue_to_surface_criterion(p1, p2))
      return;

  // Ignore virtual particles
  if (p1.is_virtual() or p2.is_virtual())
    return;

  // Check, if there's already a bond between the particles
  if (pair_bond_exists_on(p1.bonds(), p2.id(), collision_params.bond_centers))
    return;

  if (pair_bond_exists_on(p2.bonds(), p1.id(), collision_params.bond_centers))
    return;

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  // do not create bond between ghost particles
  if (p1.is_ghost() and p2.is_ghost()) {
    return;
  }
  queue_collision(p1.id(), p2.id());
}

#endif // COLLISION_DETECTION
double collision_detection_cutoff();
double collision_detection_cutoff(std::vector<int> types);
