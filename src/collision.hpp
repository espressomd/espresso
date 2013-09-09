/*
  Copyright (C) 2011,2012,2013 The ESPResSo project
  
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

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "virtual_sites_relative.hpp"
#include "virtual_sites.hpp"
#include "integrate.hpp"

#ifdef COLLISION_DETECTION

/** \name bits of possible modes for collision handling.
    To be used with \ref collision_detection_set_params.
    The modes can be combined by or-ing together. Not all combinations are possible.
    COLLISION_MODE_ERROR|COLLISION_MODE_BOND.
*/
/*@{*/
/** raise a background error on collision, to allow further processing in Tcl.
    Can be combined with a bonding mode, if desired
 */
#define COLLISION_MODE_EXCEPTION 1
/// just create bond between centers of colliding particles
#define COLLISION_MODE_BOND  2
/** create a bond between the centers of the colloiding particles,
    plus two virtual sites at the point of collision and bind them
    together. This prevents the particles from sliding against each
    other. Requires VIRTUAL_SITES_RELATIVE and COLLISION_MODE_BOND*/
#define COLLISION_MODE_VS    4
/*@}*/

typedef struct {
  /// bond type used between centers of colliding particles
  int bond_centers;
  /// bond type used between virtual sites 
  int bond_vs;
  /// particle type for virtual sites created on collision
  int vs_particle_type;
  /// collision handling mode, a combination of constants COLLISION_MODE_*
  int mode;
  /// distance at which particles are bound
  double distance;
} Collision_parameters;

/// Parameters for collision detection
extern Collision_parameters collision_params;

/** Detect a collision between two particles. In case of collision,
    a bond between the particles is added as marker and the collision is
    recorded in the queue for later processing.
*/
void detect_collision(Particle* p1, Particle* p2);

void prepare_collision_queue();

/// Handle the collisions recorded in the queue
void handle_collisions();

/** set the parameters for the collision detection
    @param mode is a bitset out of the COLLISION_MODE_* bits
    @param d is the collision distance, below that a bond is generated
    @param bond_centers is the type of the bond between the real particles
    @param bond_vs is the type of the bond between the virtual particles,
    if using noslip bonds
    @param t is the type of the virtual sites, if using noslip bonds
 */
int collision_detection_set_params(int mode, double d, int bond_centers, int bond_vs,int t);

#endif

#endif
