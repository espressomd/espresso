/*
  Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _COLLISION_H
#define _COLLISION_H

#define COLLISION_MODE_OFF 0
/// just create bond between centers of colliding particles
#define COLLISION_MODE_BOND 2
/** create a bond between the centers of the colloiding particles,
    plus two virtual sites at the point of collision and bind them
    together. This prevents the particles from sliding against each
    other. Requires VIRTUAL_SITES_RELATIVE and COLLISION_MODE_BOND*/
#define COLLISION_MODE_VS 4
/** Glue a particle to a speciffic spot on the surface of an other */
#define COLLISION_MODE_GLUE_TO_SURF 8
/// Three particle binding mode
#define COLLISION_MODE_BIND_THREE_PARTICLES 16
/*@}*/

class Collision_parameters {
public:
  Collision_parameters()
      : mode(COLLISION_MODE_OFF), distance(0.), bond_centers(-1), bond_vs(-1),
        bond_three_particles(-1){};

  /// collision handling mode, a combination of constants COLLISION_MODE_*
  int mode;
  /// distance at which particles are bound
  double distance;

  /// bond type used between centers of colliding particles
  int bond_centers;
  /// bond type used between virtual sites
  int bond_vs;
  /// particle type for virtual sites created on collision
  int vs_particle_type;

  /** Raise exception on collision */
  bool exception_on_collision;

  /// For mode "glue to surface": The distance from the particle which is to be
  /// glued to the new virtual site
  double dist_glued_part_to_vs;
  /// For mode "glue to surface": The particle type being glued
  int part_type_to_be_glued;
  /// For mode "glue to surface": The particle type to which the virtual site is
  /// attached
  int part_type_to_attach_vs_to;
  /// Particle type to which the newly glued particle is converd
  int part_type_after_glueing;
  /// first bond type (for zero degrees)) used for the three-particle bond
  /// (angle potential)
  int bond_three_particles;
  /// Number of angle bonds to use (angular resolution)
  /// different angle bonds with different equilibrium angles
  /// Are expected to have ids immediately following to bond_three_particles
  int three_particle_angle_resolution;
  /** Placement of virtual sites for MODE_VS. 0=on same particle as related to,
   * 1=on collision partner. 0.5=in the middle between */
  double vs_placement;
};
/// Parameters for collision detection
extern Collision_parameters collision_params;

#include "integrate.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "virtual_sites.hpp"

#ifdef COLLISION_DETECTION

/** \name bits of possible modes for collision handling.
    To be used with \ref collision_detection_set_params.
    The modes can be combined by or-ing together. Not all combinations are
   possible.
    COLLISION_MODE_ERROR|COLLISION_MODE_BOND.
*/
/*@{*/
/** raise a background error on collision, to allow further processing in Tcl.
    Can be combined with a bonding mode, if desired
 */

void prepare_local_collision_queue();

/// Handle the collisions recorded in the queue
void handle_collisions();

/** @brief Validates collision parameters and creates particle types if needed
 */
bool validate_collision_parameters();

/** @brief add the collision between the given particle ids to the collision
 * queue */
void queue_collision(const int part1, const int part2);

/** @brief Check additional criteria for the glue_to_surface collision mode */
inline bool glue_to_surface_criterion(const Particle *const p1,
                                      const Particle *const p2) {
  return (((p1->p.type == collision_params.part_type_to_be_glued) &&
           (p2->p.type == collision_params.part_type_to_attach_vs_to)) ||
          ((p2->p.type == collision_params.part_type_to_be_glued) &&
           (p1->p.type == collision_params.part_type_to_attach_vs_to)));
}

/** @brief Detect (and queue) a collision between the given particles. */
inline void detect_collision(const Particle *const p1, const Particle *const p2,
                             const double &dist_betw_part) {
  if (dist_betw_part > collision_params.distance)
    return;

  // If we are in the glue to surface mode, check that the particles
  // are of the right type
  if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)
    if (!glue_to_surface_criterion(p1, p2))
      return;

#ifdef VIRTUAL_SITES_RELATIVE
  // Ignore virtual particles
  if ((p1->p.isVirtual) || (p2->p.isVirtual))
    return;
#endif

  // Check, if there's already a bond between the particles
  if (pair_bond_exists_on(p1, p2, collision_params.bond_centers))
    return;

  if (pair_bond_exists_on(p2, p1, collision_params.bond_centers))
    return;

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  // do not create bond between ghost particles
  if (p1->l.ghost && p2->l.ghost) {
    return;
  }
  queue_collision(p1->p.identity, p2->p.identity);
}

#endif

inline double collision_detection_cutoff() {
#ifdef COLLISION_DETECTION
  if (collision_params.mode)
    return collision_params.distance;
#endif
  return 0.;
}

#endif
