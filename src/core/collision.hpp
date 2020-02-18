/*
 * Copyright (C) 2011-2019 The ESPResSo project
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
#ifndef _COLLISION_H
#define _COLLISION_H

/** \name bits of possible modes for collision handling.
 *  The modes can be combined by or-ing together, but not all combinations are
 *  possible.
 */
/*@{*/
#define COLLISION_MODE_OFF 0
/// Just create bond between centers of colliding particles
#define COLLISION_MODE_BOND 2
/** Create a bond between the centers of the colliding particles,
 *  plus two virtual sites at the point of collision and bind them
 *  together. This prevents the particles from sliding against each
 *  other. Requires VIRTUAL_SITES_RELATIVE and @ref COLLISION_MODE_BOND
 */
#define COLLISION_MODE_VS 4
/** Glue a particle to a specific spot on the surface of an other */
#define COLLISION_MODE_GLUE_TO_SURF 8
/// Three particle binding mode
#define COLLISION_MODE_BIND_THREE_PARTICLES 16
/*@}*/

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "integrate.hpp"
#include "random.hpp"
#include "virtual_sites.hpp"

class Collision_parameters {
public:
  Collision_parameters()
      : active(0), distance(0.), distance2(0.), rate(0.),
        bond_type(-1), vs_bond_type(-1), bond_three_particles(-1){};

  /// collision handling mode, a combination of constants COLLISION_MODE_*
  int active;

  /** Raise exception on collision */
  bool exception_on_collision;

  /// First bond type (for zero degrees) used for the three-particle bond
  /// (angle potential)
  int bond_three_particles;
  /// Number of angle bonds to use (angular resolution)
  /// different angle bonds with different equilibrium angles
  /// Are expected to have ids immediately following to bond_three_particles
  int three_particle_angle_resolution;

  double distance;

  double distance2;

  double rate;

  std::vector<int> particle_type;

  std::vector<int> particle_type_after_collision;

  std::vector<int> vs_particle_type;

  std::vector<double> distance_vs_particle;

  int bond_type;

  int vs_bond_type;
};

/// Parameters for collision detection
extern Collision_parameters collision_params;

#ifdef COLLISION_DETECTION

void prepare_local_collision_queue();

/// Handle the collisions recorded in the queue
void handle_collisions();

/** @brief Validates collision parameters and creates particle types if needed
 */
bool validate_collision_parameters();

/** @brief Add the collision between the given particle ids to the collision
 *  queue
 */
void queue_collision(std::vector<int> particles);

inline bool particle_type_criterion(Particle const &p1, Particle const &p2) {
  int count_multiples = 1;
  if(p1.p.type == p2.p.type)
    count_multiples = 2;
  if(std::count(collision_params.particle_type.begin(),
                collision_params.particle_type.end(),p1.p.type) == count_multiples)
    return (count_multiples == 2 || std::count(collision_params.particle_type.begin(),
                                               collision_params.particle_type.end(),p2.p.type) > 0);
  return 0;
}

inline bool collision_detection_criterion(Particle const &p1, Particle const &p2){
  if(collision_params.rate > 0)
    return d_random() < collision_params.rate * time_step;
  /* TODO implement other criteria */
}

inline bool virtual_site_criterion(Particle const &p1, Particle const &p2){
  return ((p1.p.is_virtual or p2.p.is_virtual) and !collision_params.vs_particle_type.empty());
}

/** @brief Detect (and queue) a collision between the given particles. */
inline void detect_collision(Particle const &p1, Particle const &p2,
                             const double &dist_betw_part2) {
  if (dist_betw_part2 > collision_params.distance2)
    return;

  // Check, if the particle types match the criteria
  if (!particle_type_criterion(p1,p2))
    return;

  // Check, if there's already a bond between the particles
  if (pair_bond_exists_on(p1, p2, collision_params.bond_type))
    return;

  if (pair_bond_exists_on(p2, p1, collision_params.bond_type))
    return;

  if(!collision_detection_criterion(p1, p2))
    return;

  if(virtual_site_criterion(p1,p2))
    return;

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  // do not create bond between ghost particles
  if (p1.l.ghost && p2.l.ghost) {
    return;
  }
  queue_collision({p1.p.identity, p2.p.identity});
}

#endif

inline double collision_detection_cutoff() {
#ifdef COLLISION_DETECTION
  return collision_params.distance2;
#endif
  return 0.;
}

#endif
