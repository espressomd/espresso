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

#include "collision.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"

// #define DEBUG

#ifdef DEBUG
#define TRACE(a) a
#else
#define TRACE(a)
#endif

#ifdef COLLISION_DETECTION

/// Data type holding the info about a single collision
typedef struct {
  int pp1; // 1st particle id
  int pp2; // 2nd particle id
  double point_of_collision[3]; 
} collision_struct;

// During force calculation, colliding particles are recorded in thequeue
// The queue is processed after force calculation, when it is save to add
// particles
static collision_struct *collision_queue;
// Number of collisions recoreded in the queue
static int number_of_collisions;

/// Parameters for collision detection
Collision_parameters collision_params = { 0, };

int collision_detection_set_params(int mode, double d, int bond_centers,
				   int bond_vs, int t)
{
  if (mode & COLLISION_MODE_VS)
    mode |= COLLISION_MODE_BOND;

  // If we don't have virtual sites, virtual site binding isn't possible.
#ifndef VIRTUAL_SITES_RELATIVE
  if (mode & COLLISION_MODE_VS)
    return 1;
#endif

  // Binding so far only works on a single cpu
  if ((mode & COLLISION_MODE_VS) && n_nodes != 1)
    return 2;

  // Check if bonded ia exist
  if ((mode & COLLISION_MODE_BOND) &&
      (bond_centers >= n_bonded_ia))
    return 3;
  if ((mode & COLLISION_MODE_VS) &&
      (bond_vs >= n_bonded_ia))
    return 3;
  
  if ((mode & COLLISION_MODE_BOND) &&
      (bonded_ia_params[bond_centers].num != 1))
    return 4;
  
  if ((mode & COLLISION_MODE_VS) && !(bonded_ia_params[bond_vs].num == 1 ||
				      bonded_ia_params[bond_vs].num == 2))
    return 5;

  // Set params
  collision_params.mode=mode;
  collision_params.bond_centers=bond_centers;
  collision_params.bond_vs=bond_vs;
  collision_params.distance=d;
  collision_params.vs_particle_type=t;

  make_particle_type_exist(t);

  mpi_bcast_collision_params();

  return 0;
}

// Detect a collision between the given particles.
// Add it to the queue in case virtual sites should be added at the point of collision
void detect_collision(Particle* p1, Particle* p2)
{
  // The check, whether collision detection is actually turned on is performed in forces.hpp

  double dist_betw_part, vec21[3]; 
  int part1, part2, size;

  TRACE(printf("%d: consider particles %d and %d\n", this_node, p1->p.identity, p2->p.identity));

  // Obtain distance between particles
  dist_betw_part = sqrt(distance2vec(p1->r.p, p2->r.p, vec21));
  if (dist_betw_part > collision_params.distance)
    return;

  TRACE(printf("%d: particles %d and %d on bonding distance %lf\n", this_node, p1->p.identity, p2->p.identity, dist_betw_part));

  part1 = p1->p.identity;
  part2 = p2->p.identity;
      
  // Retrieving the particles from local_particles is necessary, because the particle might be a
  // ghost, and those can't store bonding info.
  p1 = local_particles[part1];
  p2 = local_particles[part2];

#ifdef VIRTUAL_SITES_RELATIVE
  // Ignore virtual particles
  if ((p1->p.isVirtual) || (p2->p.isVirtual))
    return;
#endif

  // Check, if there's already a bond between the particles
  // First check the bonds of p1
  if (p1->bl.e) {
    int i = 0;
    while(i < p1->bl.n) {
      size = bonded_ia_params[p1->bl.e[i]].num;
      
      if (p1->bl.e[i] == collision_params.bond_centers &&
          p1->bl.e[i + 1] == part2) {
        // There's a bond, already. Nothing to do for these particles
        return;
      }
      i += size + 1;
    }
  }
  if (p2->bl.e) {
    // Check, if a bond is already stored in p2
    int i = 0;
    while(i < p2->bl.n) {
      size = bonded_ia_params[p2->bl.e[i]].num;

      /* COMPARE P2 WITH P1'S BONDED PARTICLES*/

      if (p2->bl.e[i] == collision_params.bond_centers &&
          p2->bl.e[i + 1] == part1) {
        return;
      }
      i += size + 1;
    }
  }

  TRACE(printf("%d: no previous bond, binding\n", this_node));

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  /* create marking bond between the colliding particles immediately */
  if (collision_params.mode & COLLISION_MODE_BOND) {
    int bondG[2];
    int primary = part1, secondary = part2;
    // put the bond to the physical particle; at least one partner always is
    if (p1->l.ghost) {
      primary = part2;
      secondary = part1;
    }
    bondG[0]=collision_params.bond_centers;
    bondG[1]=secondary;
    local_change_bond(primary, bondG, 0);
  }

  if (collision_params.mode & (COLLISION_MODE_VS | COLLISION_MODE_EXCEPTION)) {
    /* If we also create virtual sites or throw an exception, we add the collision
       to the queue to process later */

    // Point of collision
    double new_position[3];
    for (int i=0;i<3;i++) {
      new_position[i] = p1->r.p[i] - vec21[i] * 0.50;
    }
       
    number_of_collisions++;

    // Allocate mem for the new collision info
    collision_queue = (collision_struct *) realloc (collision_queue, (number_of_collisions) * sizeof(collision_struct));
      
    // Save the collision      
    collision_queue[number_of_collisions-1].pp1 = part1;
    collision_queue[number_of_collisions-1].pp2 = part2;
    for (int i=0;i<3;i++) {
      collision_queue[number_of_collisions-1].point_of_collision[i] = new_position[i]; 
    }
  }
}

void prepare_collision_queue()
{
  
  number_of_collisions=0;

  collision_queue = (collision_struct *) malloc (sizeof(collision_struct));

}

// Handle the collisions stored in the queue
void handle_collisions ()
{
  // printf("number of collisions in handle collision are %d\n",number_of_collisions);

  for (int i = 0; i < number_of_collisions; i++) {
    //  printf("Handling collision of particles %d %d\n", collision_queue[i].pp1, collision_queue[i].pp2);
    //  fflush(stdout);

    if (collision_params.mode & (COLLISION_MODE_EXCEPTION)) {

      // if desired, raise a runtime exception (background error) on collision
      char *exceptiontxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
      int id1, id2;
      if (collision_queue[i].pp1 > collision_queue[i].pp2) {
	id1 = collision_queue[i].pp2;
	id2 = collision_queue[i].pp1;
      }
      else {
	id1 = collision_queue[i].pp1;
	id2 = collision_queue[i].pp2;
      }
      ERROR_SPRINTF(exceptiontxt, "{collision between particles %d and %d} ",
		    id1, id2);
    }

    /* If we don't have virtual_sites_relative, only bonds between centers of 
       colliding particles are possible and nothing is to be done here */
#ifdef VIRTUAL_SITES_RELATIVE
    if (collision_params.mode & COLLISION_MODE_VS) {

      // add virtual sites placed at the point of collision and bind them
      int bondG[3];
  
      // Virtual site related to first particle in the collision
      place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
      vs_relate_to(max_seen_particle,collision_queue[i].pp1);
      (local_particles[max_seen_particle])->p.isVirtual=1;
      (local_particles[max_seen_particle])->p.type=collision_params.vs_particle_type;
  
      // Virtual particle related to 2nd particle of the collision
      place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
      vs_relate_to(max_seen_particle,collision_queue[i].pp2);
      (local_particles[max_seen_particle])->p.isVirtual=1;
      (local_particles[max_seen_particle])->p.type=collision_params.vs_particle_type;
  
      switch (bonded_ia_params[collision_params.bond_vs].num) {
      case 1: {
	// Create bond between the virtual particles
	bondG[0] = collision_params.bond_vs;
	bondG[1] = max_seen_particle-1;
	local_change_bond(max_seen_particle, bondG, 0);
	break;
      }
      case 2: {
	// Create 1st bond between the virtual particles
	bondG[0] = collision_params.bond_vs;
	bondG[1] = collision_queue[i].pp1;
	bondG[2] = collision_queue[i].pp2;
	local_change_bond(max_seen_particle,   bondG, 0);
	local_change_bond(max_seen_particle-1, bondG, 0);
	break;
      }
      }
    }
#endif
  }

  // Reset the collision queue
  number_of_collisions = 0;
  free(collision_queue);

  announce_resort_particles();
}

#endif
