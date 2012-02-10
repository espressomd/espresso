/*
  Copyright (C) 2011,2012 The ESPResSo project
  
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

#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "virtual_sites.h"
#include "integrate.h"

#ifdef COLLISION_DETECTION

/// Data type holding the info about a single collision
typedef struct {
  int pp1; // 1st particle id
  int pp2; // 2nd particle id
  double point_of_collision[3]; 
} collision_struct;

/// Bond type used between centers of colliding particles
extern int collision_detection_bond_centers;

/// bond type used between virtual sites 
extern int collision_detection_bond_vs;

/// Particle type for virtual sites created on collision
extern int collision_vs_particle_type;

/** Bonding mode: 
    0 = off
    1 = only create bond between centers of colliiding particles
    2 = also create two virtual sites at the point o collision and bind the 
    together. This prevents the particles from sliding against each other. 
    (Requires VIRTUAL_SITES_RELATIVE) */
extern int collision_detection_mode;

/// Distance at which particles are bound
extern double collision_distance;

/** Collision detection mod
    0=off
    1=bind centers
    2=bind at point of collision
    Never write to this variable. Use collision_detection_set_params()
*/
extern int collision_detection_mode;

/** Detect a collision between two particles. In case of collision,
    a bond between the particles is added as marker and the collision is
    recorded in the queue for later processing.
*/
void detect_collision(Particle* p1, Particle* p2);

void prepare_collision_queue();

/// Handle the collisions recorded in the queue
void handle_collisions();

/** set the parameters for the collision detection
    @param mode is 0, 1 or 2 for off, simple bond or noslip
    @param d is the collision distance, below that a bond is generated
    @param bond_centers is the type of the bond between the real particles
    @param bond_vs is the type of the bond between the virtual particles,
    if using noslip bonds
    @param t is the type of the virtual sites, if using noslip bonds
 */
int collision_detection_set_params(int mode, double d, int bond_centers, int bond_vs,int t);

#endif

#endif
