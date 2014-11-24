/*
  Copyright (C) 2011,2012,2013,2014 The ESPResSo project
  
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
#include "grid.hpp"
#include "domain_decomposition.hpp"


using namespace std;

#define DEBUG

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
static collision_struct *gathered_queue;

// Number of collisions recoreded in the queue
static int number_of_collisions, total_collisions;

/// Parameters for collision detection
Collision_parameters collision_params = { 0, };

int collision_detection_set_params(int mode, double d, int bond_centers, int bond_vs,int t, int bond_three_particles, int angle_resolution)
{
  if (mode & COLLISION_MODE_VS)
    mode |= COLLISION_MODE_BOND;

  if (mode & COLLISION_MODE_BIND_THREE_PARTICLES)
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

  // Check that the bonds have the right number of partners
  if ((mode & COLLISION_MODE_BOND) &&
      (bonded_ia_params[bond_centers].num != 1))
    return 4;
  
  if ((mode & COLLISION_MODE_VS) && !(bonded_ia_params[bond_vs].num == 1 ||
				      bonded_ia_params[bond_vs].num == 2))
    return 5;
  
  if (mode & COLLISION_MODE_BIND_THREE_PARTICLES) {
    if (bond_three_particles + angle_resolution >= n_bonded_ia)
      return 6;
    
    for (int i = bond_three_particles; i <= bond_three_particles + angle_resolution; i++) {
      if (bonded_ia_params[i].num != 2)
        return 7;
    }
  }

  // Set params
  collision_params.mode=mode;
  collision_params.bond_centers=bond_centers;
  collision_params.bond_vs=bond_vs;
  collision_params.distance=d;
  collision_params.vs_particle_type=t;
  collision_params.bond_three_particles=bond_three_particles;
  collision_params.three_particle_angle_resolution=angle_resolution;

  if (mode & COLLISION_MODE_VS)
    make_particle_type_exist(t);

  mpi_bcast_collision_params();

  return 0;
}

// Detect a collision between the given particles.
// Add it to the queue in case virtual sites should be added at the point of collision
void detect_collision(Particle* p1, Particle* p2)
{
  // The check, whether collision detection is actually turned on is performed in forces.hpp

  int part1, part2, size;
  int counts[n_nodes];
  //TRACE(printf("%d: consider particles %d and %d\n", this_node, p1->p.identity, p2->p.identity));

  double vec21[3];
  // Obtain distance between particles
  double dist_betw_part = sqrt(distance2vec(p1->r.p, p2->r.p, vec21));
  if (dist_betw_part > collision_params.distance)
    return;

  //TRACE(printf("%d: particles %d and %d on bonding distance %lf\n", this_node, p1->p.identity, p2->p.identity, dist_betw_part));

  part1 = p1->p.identity;
  part2 = p2->p.identity;
      
  // Retrieving the particles from local_particles is necessary, because the particle might be a
  // ghost, and those can't store bonding info.
  p1 = local_particles[part1];
  p2 = local_particles[part2];

#ifdef COLLISION_USE_BROKEN_PARALLELIZATION
  // Ignore particles too distant to be on the same processor
  if (!p1 || !p2)
    return; 
#endif

#ifdef VIRTUAL_SITES_RELATIVE
  // Ignore virtual particles
  if ((p1->p.isVirtual) || (p2->p.isVirtual))
    return;
#endif

  if (p1==p2)
    return;

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

  //TRACE(printf("%d: no previous bond, binding\n", this_node));

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  /* create marking bond between the colliding particles immediately */
  if (collision_params.mode & COLLISION_MODE_BOND) {
    int bondG[2];
    int primary = part1, secondary = part2;

    // do not create bond between ghost particles
    if (p1->l.ghost && p2->l.ghost) {
       TRACE(printf("Both particles %d and %d are ghost particles", p1->p.identity, p2->p.identity));
       return;
    }
    // put the bond to the physical particle; at least one partner always is
    if (p1->l.ghost) {
      primary = part2;
      secondary = part1;
      //TRACE(printf("%d: particle-%d is ghost", this_node, p1->p.identity));
    }
    bondG[0]=collision_params.bond_centers;
    bondG[1]=secondary;
    local_change_bond(primary, bondG, 0);
  }

  if (collision_params.mode & (COLLISION_MODE_VS | COLLISION_MODE_EXCEPTION | COLLISION_MODE_BIND_THREE_PARTICLES)) {
    /* If we also create virtual sites or bind three particles, or throw an exception, we add the collision
       to the queue to process later */

    // Point of collision
    double new_position[3];
    for (int i=0;i<3;i++) {
      new_position[i] = p1->r.p[i] - vec21[i] * 0.50;
    }
    


    number_of_collisions++;
    //counts[this_node]=number_of_collisions;
    //TRACE(printf("%d: particles %d and %d with collision of number: %d and counts: %d\n", this_node, p1->p.identity, p2->p.identity, number_of_collisions, counts[this_node]));


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
  total_collisions=0;
  collision_queue = (collision_struct *) malloc (sizeof(collision_struct));

  gathered_queue = (collision_struct *) malloc (sizeof(collision_struct));
}

// See comments in handle_collsion_queue()
void coldet_do_three_particle_bond(Particle* p, Particle* p1, Particle* p2)
{
  double vec21[3];
  // If p1 and p2 are not closer or equal to the cutoff distance, skip
  // p1:
  get_mi_vector(vec21,p->r.p,p1->r.p);
  if (sqrt(sqrlen(vec21)) > collision_params.distance)
    return;
  // p2:
  get_mi_vector(vec21,p->r.p,p2->r.p);
  if (sqrt(sqrlen(vec21)) > collision_params.distance)
    return;

  //TRACE(printf("%d: checking three particle bond %d %d %d\n", this_node, p1->p.identity, p->p.identity, p2->p.identity));

  // Check, if there already is a three-particle bond centered on p 
  // with p1 and p2 as partners. If so, skip this triplet.
  // Note that the bond partners can appear in any order.
 
  // Iterate over existing bonds of p

  if (p->bl.e) {
    int b = 0;
    while (b < p->bl.n) {
      int size = bonded_ia_params[p->bl.e[b]].num;

      //TRACE(printf("%d:--1-- checking bond of type %d and length %d of particle %d\n", this_node, p->bl.e[b], bonded_ia_params[p->bl.e[b]].num, p->p.identity));
 
      if (size==2) {
        // Check if the bond type is within the range used by the collision detection,
        if ((p->bl.e[b] >= collision_params.bond_three_particles) & (p->bl.e[b] <=collision_params.bond_three_particles + collision_params.three_particle_angle_resolution)) {
          // check, if p1 and p2 are the bond partners, (in any order)
          // if yes, skip triplet
          if (
              ((p->bl.e[b+1]==p1->p.identity) && (p->bl.e[b+2] ==p2->p.identity))
              ||
              ((p->bl.e[b+1]==p2->p.identity) && (p->bl.e[b+2] ==p1->p.identity))
              )
            return;
        } // if bond type 
      } // if size==2
      
      // Go to next bond
      b += size + 1;
    } // bond loop
  } // if bond list defined

  TRACE(printf("%d: proceeding to install three particle bond %d %d %d\n", this_node, p1->p.identity, p->p.identity, p2->p.identity));

  // If we are still here, we need to create angular bond
  // First, find the angle between the particle p, p1 and p2
  double cosine=0.0;
  
  double vec1[3],vec2[3];
  /* vector from p to p1 */
  get_mi_vector(vec1, p->r.p, p1->r.p);
  // Normalize
  double dist2 = sqrlen(vec1);
  double d1i = 1.0 / sqrt(dist2);
  for(int j=0;j<3;j++) vec1[j] *= d1i;
  
  /* vector from p to p2 */
  get_mi_vector(vec2, p->r.p, p2->r.p);
  // normalize
  dist2 = sqrlen(vec2);
  double d2i = 1.0 / sqrt(dist2);
  for(int j=0;j<3;j++) vec2[j] *= d2i;
  
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  
  // Handle case where cosine is nearly 1 or nearly -1
  if ( cosine >  TINY_COS_VALUE)  
    cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  
    cosine = -TINY_COS_VALUE;
  
  // Bond angle
  double phi =  acos(cosine);
  
  // We find the bond id by dividing the range from 0 to pi in 
  // three_particle_angle_resolution steps and by adding the id
  // of the bond for zero degrees.
  int bond_id =floor(phi/M_PI * (collision_params.three_particle_angle_resolution-1) +0.5) +collision_params.bond_three_particles;
  
  // Create the bond
  
  // First, fill bond data structure
  int bondT[3];
  bondT[0] = bond_id;
  bondT[1] = p1->p.identity;
  bondT[2] = p2->p.identity;
  local_change_bond(p->p.identity, bondT, 0);
  
}

// Handle the collisions stored in the queue
void handle_collisions ()
{
  for (int i = 0; i < number_of_collisions; i++) {
    //TRACE(printf("Handling collision of particles %d %d\n", collision_queue[i].pp1, collision_queue[i].pp2));

    if (collision_params.mode & (COLLISION_MODE_EXCEPTION)) {

      int id1, id2;
      if (collision_queue[i].pp1 > collision_queue[i].pp2) {
	id1 = collision_queue[i].pp2;
	id2 = collision_queue[i].pp1;
      }
      else {
	id1 = collision_queue[i].pp1;
	id2 = collision_queue[i].pp2;
      }
      ostringstream msg;
      msg << "collision between particles " << id1 << " and " <<id2;
      runtimeError(msg);
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

  // three-particle-binding part


  int counts[n_nodes];                          // number of collisions on each proc
  //int total_collisions;                         // total number of collisions
  //collision_struct gathered_queue[10000];		// collisions on all nodes
  int displacements[n_nodes];                   // offsets into collisions
  
    MPI_Datatype collisiontype, oldtypes[2];
    int blockcounts[2];
    // MPI_Aint type used to be consistent with syntax of MPI_Type_extent routine
    MPI_Aint offsets[2], extent;
    
    // Setup description of the 2 MPI_INT fields pp1, pp2
    offsets[0] = 0;
    oldtypes[0] = MPI_INT;
    blockcounts[0] = 2;
    
    // Setup description of the 1 MPI_DOUBLE field position
    // Need to first figure offset by getting size of MPI_INT
    MPI_Type_extent(MPI_INT, &extent);
    offsets[1] = 2*extent;
    oldtypes[1] = MPI_DOUBLE;
    blockcounts[1] = 3;
    
    // Now define structured type and commit it
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &collisiontype);
    MPI_Type_commit(&collisiontype);

  // Initialize number of collisions gathered from all processors
  for (int a=0;a<n_nodes;a++)
      counts[a]=0;
    
  // Total number of collisions
  MPI_Allreduce(&number_of_collisions, &total_collisions, 1, MPI_INT, MPI_SUM, comm_cart);

  // Gather number of collisions
  MPI_Allgather(&number_of_collisions, 1, MPI_INT, counts, 1, MPI_INT, comm_cart);

  // initialize displacement information for all nodes
  displacements[0]=0;
  
  // Find where to place collision information for each processor
  for (int k=1; k<n_nodes; k++)
      displacements[k]=displacements[k-1]+counts[k-1];
    
  //TRACE(printf("counts [%d] = %d and number of collisions = %d and diplacements = %d and total collisions = %d\n", this_node, counts[this_node], number_of_collisions, displacements[this_node], total_collisions));
    
    // Allocate mem for the new collision info
    gathered_queue = (collision_struct *) realloc (gathered_queue, (total_collisions) * sizeof(collision_struct));

       MPI_Allgatherv(collision_queue, counts[this_node], collisiontype, gathered_queue, counts, displacements, collisiontype, comm_cart);

int part1, part2;
for (int k=0;k<total_collisions;k++) {
part1=gathered_queue[k].pp1;
part2=gathered_queue[k].pp2;
//TRACE(printf("%d: collisions of %d %d\n", this_node, part1, part2));
}

  //if (number_of_collisions>0) {

   if (collision_params.mode & (COLLISION_MODE_BIND_THREE_PARTICLES)) {

   // If we don't have domain decomposition, we need to do a full sweep over all
   // particles in the system. (slow)
   if (cell_structure.type!=CELL_STRUCTURE_DOMDEC)
   {
       Cell *cell;
       Particle *p, *p1, *p2;
       // first iterate over cells, get one of the cells and find how many particles in this cell
       for (int c=0; c<local_cells.n; c++) {
           cell=local_cells.cell[c];
           // iterate over particles in the cell
           for (int a=0; a<cell->n; a++) {
               p=&cell->part[a];
               // for all p:
               for (int ij=0; ij<total_collisions; ij++) {
                   p1=local_particles[gathered_queue[ij].pp1];
                   p2=local_particles[gathered_queue[ij].pp2];
  
  		 // Check, whether p is equal to one of the particles in the
  		 // collision. If so, skip
  		 if ((p->p.identity ==p1->p.identity) || ( p->p.identity == p2->p.identity)) {
  		   continue;
  		 }
  
                   // The following checks, 
  		 // if the particle p is closer that the cutoff from p1 and/or p2.
  		 // If yes, three particle bonds are created on all particles
  		 // which have two other particles within the cutoff distance,
  		 // unless such a bond already exists
  		 
  		 // We need all cyclical permutations, here 
  		 // (bond is placed on 1st particle, order of bond partners
  		 // does not matter, so we don't neet non-cyclic permutations):
                   coldet_do_three_particle_bond(p,p1,p2);
                   coldet_do_three_particle_bond(p1,p,p2);
                   coldet_do_three_particle_bond(p2,p,p1);
  
               }
           }
       }
  
    } // if cell structure = domain decomposition
    else
    { // We have domain decomposition
    
    // Indices of the cells in which the colliding particles reside
    int cellIdx[2][3];
    
    // Iterate over collision queue
//if (this_node==0) {
    for (int id=0;id<total_collisions;id++) {

   //TRACE(printf("domain decomposition\n"));

      // Get first cell Idx
      Particle* p1=local_particles[gathered_queue[id].pp1];
      Particle* p2=local_particles[gathered_queue[id].pp2];
      dd_position_to_cell_indices(p1->r.p,cellIdx[0]);
      dd_position_to_cell_indices(p2->r.p,cellIdx[1]);

//TRACE(printf("%d: We have particles from gathered_queue %d %d and found first cell idx\n", this_node, p1->p.identity, p2->p.identity));

      // Iterate over the cells + their neighbors
      // if p1 and p2 are in the same cell, we don't need to consider it 2x
      int lim=1;

      if ((cellIdx[0][0]==cellIdx[1][0]) && (cellIdx[0][1]==cellIdx[1][1]) && (cellIdx[0][2]==cellIdx[1][2]))
        lim=0; // Only consider the 1st cell
//TRACE(printf("%d:Step 1 wit particles %d %d\n", this_node, p1->p.identity, p2->p.identity));
      for (int j=0;j<=lim;j++)
      {
       // Iterate the cell with indices cellIdx[j][] and all its neighbors.
       // code taken from dd_init_cell_interactions()
       for(int p=cellIdx[j][0]-1; p<=cellIdx[j][0]+1; p++)	
         for(int q=cellIdx[j][1]-1; q<=cellIdx[j][1]+1; q++)
	   for(int r=cellIdx[j][2]-1; r<=cellIdx[j][2]+1; r++) {   
	    int ind2 = get_linear_index(p,q,r,dd.ghost_cell_grid);
	    Cell* cell=cells+ind2;
//TRACE(printf("%d:Step 2 wit particles %d %d\n", this_node, p1->p.identity, p2->p.identity));
	    // Iterate over particles in this cell
            for (int a=0; a<cell->n; a++) {
               Particle* P=&cell->part[a];
//TRACE(printf("%d:Step 3 with particles %d %d and %d\n", this_node, p1->p.identity, p2->p.identity, P->p.identity));
               // for all p:
  	       // Check, whether p is equal to one of the particles in the
  	       // collision. If so, skip
  	       if ((P->p.identity ==p1->p.identity) || (P->p.identity == p2->p.identity) || (p1->p.identity == p2->p.identity)) {
               //TRACE(printf("same particle\n"));
  		   continue;
  	       }
TRACE(printf("%d: I found 3th match %d to %d %d\n", this_node, P->p.identity, p1->p.identity, p2->p.identity));
               // The following checks, 
               // if the particle p is closer that the cutoff from p1 and/or p2.
               // If yes, three particle bonds are created on all particles
  	       // which have two other particles within the cutoff distance,
               // unless such a bond already exists
 //TRACE(printf("%d: Handling collision of particles %d %d %d with total collisions %d\n", this_node, p1->p.identity, P->p.identity, p2->p.identity, total_collisions));
  	       // We need all cyclical permutations, here 
               // (bond is placed on 1st particle, order of bond partners
  	       // does not matter, so we don't need non-cyclic permutations):

//if (this_node==0) {
if (P->l.ghost) {
continue;
}
 TRACE(printf("%d: Handling collision of particles FIRST CONFIGURATION %d %d %d\n", this_node, p1->p.identity, P->p.identity, p2->p.identity));

               coldet_do_three_particle_bond(P,p1,p2);

if (p1->l.ghost) {
continue;
}
 TRACE(printf("%d: Handling collision of particles SECOND CONFIGURATION %d %d %d\n", this_node, P->p.identity, p1->p.identity, p2->p.identity));

               coldet_do_three_particle_bond(p1,P,p2);

if (p2->l.ghost) {
continue;
}
 TRACE(printf("%d: Handling collision of particles THIRD CONFIGURATION %d %d %d\n", this_node, P->p.identity, p2->p.identity, p1->p.identity));

  	       coldet_do_three_particle_bond(p2,P,p1);
//} // only node-0 creates bonds
             } // loop over particles
           } // Cell loop

	 } // Loop over 1st and 2nd particle, in case they are in diferent cells

        } // Loop over collision queue
    
      } // If domain decomposition
    } // if three particle binding
  //} // if number of collisions >0 (three particle binding)

  // Reset the collision queue
  number_of_collisions = 0;
  total_collisions = 0;
  free(collision_queue);
  free(gathered_queue);

  announce_resort_particles();
}

#endif
