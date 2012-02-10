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

#include "collision.h"
#include "cells.h"
#include "communication.h" 

#ifdef COLLISION_DETECTION

// During force calculation, colliding particles are recorded in thequeue
// The queue is processed after force calculation, when it is save to add
// particles
collision_struct * collision_queue;

// Number of collisions recoreded in the queue
int number_of_collisions;

// Bond type used between centers of colliding particles
int collision_detection_bond_centers=0;

// bond type used between virtual sites 
int collision_detection_bond_vs=1;

// Particle type for virtual sites created on collision
int collision_vs_particle_type;

// Bonding mode: 
// 0 = off
// 1 = only create bond between centers of colliiding particles
// 2 = also create two virtual sites at the point o collision and bind the 
// together. This prevents the particles from sliding against each other. 
// (Requires VIRTUAL_SITES_RELATIVE)
int collision_detection_mode=0;

// Distance at which particles are bound
double collision_distance;

int collision_detection_set_params(int mode, double d, int bond_centers, int bond_vs,int t)
{
  // If we don't have virtual sites relative feature, mode 2 isn't allowed.
#ifndef VIRTUAL_SITES_RELATIVE
  if (mode==2)
    return 1;
#endif

  // Collision detection only works on a single cpu
  if (mode!=0 && n_nodes !=1)
    return 2;

  // Check if bonded ia exist
  if (((mode>0) && (bond_centers+1>n_bonded_ia)) || ((mode==2) && (bond_vs+1>n_bonded_ia))) 
    return 3;

  if (bonded_ia_params[bond_centers].num != 1)
    return 4;

  if ((mode ==2) && (
	bonded_ia_params[bond_vs].num != 1 &&
	bonded_ia_params[bond_vs].num != 2))
    return 5;

  // Set params
  collision_detection_mode=mode;
  collision_detection_bond_centers=bond_centers;
  collision_detection_bond_vs=bond_vs;
  collision_distance=d;
  collision_vs_particle_type=t;
  make_particle_type_exist(t);
  recalc_forces=1;

  return 0;
}

// Detect a collision between the given particles.
// Add it to the queue in case virtual sites should be added at the point of collision
void detect_collision(Particle* p1, Particle* p2)
{
  // The check, whether collision detection is actually turned on is performed in forces.h

  double dist_betw_part, vec21[3]; 
  int part1, part2, size;

  // Obtain distance between particles
  dist_betw_part = distance2vec(p1->r.p, p2->r.p, vec21);
  if (dist_betw_part > collision_distance)
    return;

  part1 = p1->p.identity;
  part2 = p2->p.identity;
      
  // Retrieving the particles from local_particles is necessary, because the particle might be a
  // ghost, and those don't contain bonding info
  p1 = local_particles[part1];
  p2 = local_particles[part2];

#ifdef VIRTUAL_SITES_RELATIVE
  // Ignore virtual particles
  if ((p1->p.isVirtual) || (p2->p.isVirtual))
    return;
#endif

  // Check, if there's already a bond between the particles
  // First check the bonds of p1 
  int i = 0;
  while(i < p1->bl.n) {
    size = bonded_ia_params[p1->bl.e[i]].num;

    if (p1->bl.e[i] == collision_detection_bond_centers &&
        p1->bl.e[i + 1] == part2) {
      // There's a bond, already. Nothing to do for these particles
      return;
    }
    i += size + 1;
  }
      
  // Check, if a bond is already stored in p2
  i = 0;
  while(i < p2->bl.n) {
    size = bonded_ia_params[p2->bl.e[i]].num;

    /* COMPARE P2 WITH P1'S BONDED PARTICLES*/

    if (p2->bl.e[i] == collision_detection_bond_centers &&
        p2->bl.e[i + 1] == part1) {
      return;
    }
    i += size + 1;
  }


  // If we're still here, there is no previous bond between the particles

  // Create the bond between the particles
  int bondG[2];
  bondG[0]=collision_detection_bond_centers;
  bondG[1]=part2;
  local_change_bond(part1, bondG, 0);
  
  // If we also create virtual sites, we add the collision to the que to later add vs
  if (collision_detection_mode==2)
    {
      // Insert collision info into the queue
      
      // Point of collision
      double new_position[3];
      for (i=0;i<3;i++) {
	new_position[i] =p1->r.p[i] - vec21[i] * 0.50;
      }
       
      number_of_collisions = number_of_collisions+1;
      // Allocate mem for the new collision info
      collision_queue = (collision_struct *) realloc (collision_queue, (number_of_collisions) * sizeof(collision_struct));
      
      // Save the collision      
      collision_queue[number_of_collisions-1].pp1 = part1;
      collision_queue[number_of_collisions-1].pp2 = part2;
      for (i=0;i<3;i++) {
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
  // If we don't have virtual_sites_relative, only bonds between centers of 
  // colliding particles are possible and nothing is to be done here
#ifdef VIRTUAL_SITES_RELATIVE
  // Does the user want bonds between virtual sites placed at the point of collision
  if (collision_detection_mode==2)
    {

      //	printf("number of collisions in handle collision are %d\n",number_of_collisions);  
      int bondG[3], i;

      if (number_of_collisions > 0) {
	// Go through the queue
	for (i=0;i<number_of_collisions;i++) {
	  //  printf("Handling collision of particles %d %d\n", collision_queue[i].pp1, collision_queue[i].pp2);
	  //  fflush(stdout);
   

	  // The following lines will remove the relative velocity from
	  // colliding particles
	  //   double v[3];
	  //   for (j=0;j<3;j++)
	  //   {
	  //    v[j] =0.5 *((local_particles[collision_queue[i].pp1])->m.v[j] +(local_particles[collision_queue[i].pp2])->m.v[j]);
	  //    (local_particles[collision_queue[i].pp1])->m.v[j] =v[j];
	  //    (local_particles[collision_queue[i].pp2])->m.v[j] =v[j];
	  //   }

	  // Create virtual sites and bind them together
  
	  // Virtual site related to first particle in the collision
	  place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
	  vs_relate_to(max_seen_particle,collision_queue[i].pp1);
	  (local_particles[max_seen_particle])->p.isVirtual=1;
	  (local_particles[max_seen_particle])->p.type=collision_vs_particle_type;
  
	  // Virtual particle related to 2nd particle of the collision
	  place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
	  vs_relate_to(max_seen_particle,collision_queue[i].pp2);
	  (local_particles[max_seen_particle])->p.isVirtual=1;
	  (local_particles[max_seen_particle])->p.type=collision_vs_particle_type;
  
          switch (bonded_ia_params[collision_detection_bond_vs].num) {
	  case 1: {
	    // Create bond between the virtual particles
	    bondG[0] = collision_detection_bond_vs;
	    bondG[1] = max_seen_particle-1;
	    local_change_bond(max_seen_particle, bondG, 0);
            break;
          }
	  case 2: {
	    // Create 1st bond between the virtual particles
	    bondG[0] = collision_detection_bond_vs;
	    bondG[1] = collision_queue[i].pp1;
	    bondG[2] = collision_queue[i].pp2;
	    local_change_bond(max_seen_particle,   bondG, 0);
	    local_change_bond(max_seen_particle-1, bondG, 0);
	    break;
          }
          }
	}

      }

    }
#endif

  // Reset the collision queue	 
  number_of_collisions = 0;
  free(collision_queue);

}


#endif
