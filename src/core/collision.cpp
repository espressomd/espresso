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

#include "collision.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "interaction_data.hpp" 




#ifdef COLLISION_DETECTION_DEBUG
#define TRACE(a) a
#else
#define TRACE(a)
#endif

#ifdef COLLISION_DETECTION

/// Data type holding the info about a single collision
typedef struct {
  int pp1; // 1st particle id
  int pp2; // 2nd particle id
} collision_struct;

// During force calculation, colliding particles are recorded in thequeue
// The queue is processed after force calculation, when it is save to add
// particles
static collision_struct *collision_queue = 0;
static collision_struct *gathered_queue = 0;

// Number of collisions recoreded in the queue
static int number_of_collisions, total_collisions;

/// Parameters for collision detection
Collision_parameters collision_params = { 0, };


/** @brief Return true if a bond between the centers of the colliding particles needs to be placed. At this point, all modes need this */
inline bool bind_centers() {
  return  collision_params.mode !=COLLISION_MODE_OFF;
}


bool validate_collision_parameters()
{
  // If mode is OFF, no further checks
  if (collision_params.mode ==COLLISION_MODE_OFF) {
    return true;
  }
  // Validate distance
  if (collision_params.mode != COLLISION_MODE_OFF) {
    if (collision_params.distance<=0.) {
      runtimeErrorMsg() << "collision_detection distance must be >0";
      return false;
    }
    if (collision_params.distance >min_global_cut) {
      runtimeErrorMsg() << "The minimum global cutoff (System.min_global_cut) must be larger or equal the collision detection distance.";
    }
  }



#ifndef VIRTUAL_SITES_RELATIVE
  // The collision modes involving virutal istes also requires the creation of a bond between the colliding 
  // If we don't have virtual sites, virtual site binding isn't possible.
  if ((collision_params.mode & COLLISION_MODE_VS) || (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)) {
    runtimeErrorMsg() << "Virtual sites based collisoin modes modes require the VIRTUAL_SITES feature";
    return false;
  }
#endif


  // Check vs placement parameter
#ifdef VIRTUAL_SITES
  if (collision_params.mode & COLLISION_MODE_VS) {
    if ((collision_params.vs_placement<0) || (collision_params.vs_placement >1)) {
       runtimeErrorMsg() << "The collision detection vs_placement parameter needs to be between 0 and 1."; 
       return false;
    }
  }
#endif

  // For vs based methods, Binding so far only works on a single cpu
  if ((collision_params.mode & COLLISION_MODE_VS) ||(collision_params.mode & COLLISION_MODE_GLUE_TO_SURF))
    if (n_nodes != 1) {
      runtimeErrorMsg() << "Virtual sites based collision modes only work on a single node.";
      return false;
    }

  // Check if bonded ia exist
  if ((collision_params.mode & COLLISION_MODE_BOND) &&
      (collision_params.bond_centers >= n_bonded_ia)) {
    runtimeErrorMsg() << "The bond type to be used for binding particle centers does not exist"; 
    return false;
  }
  
  if ((collision_params.mode & COLLISION_MODE_VS) &&
      (collision_params.bond_vs >= n_bonded_ia)) {
    runtimeErrorMsg() << "The bond type to be used for binding virtual sites does not exist"; 
    return false;
  }
  
  // If the bond type to bind particle centers is not a pair bond...
  // Check that the bonds have the right number of partners
  if ((collision_params.mode & COLLISION_MODE_BOND) &&
      (bonded_ia_params[collision_params.bond_centers].num != 1)) {
    runtimeErrorMsg() << "The bond type to be used for binding particle centers needs to be a pair bond"; 
    return false;
  }
  
  // The bond between the virtual sites can be pair or triple
  if ((collision_params.mode & COLLISION_MODE_VS) && !(bonded_ia_params[collision_params.bond_vs].num == 1 ||
				      bonded_ia_params[collision_params.bond_vs].num == 2)) {
    runtimeErrorMsg() << "The bond type to be used for binding virtual sites needs to be a pair or three-particle bond"; 
    return false;
  }
  
  if (collision_params.mode & COLLISION_MODE_BIND_THREE_PARTICLES) {
    if (collision_params.bond_three_particles + collision_params.three_particle_angle_resolution > n_bonded_ia) {
      runtimeErrorMsg() << "Insufficient bonds defined for three particle binding.";
      return false;
    }

    for (int i = collision_params.bond_three_particles; i < collision_params.bond_three_particles + collision_params.three_particle_angle_resolution; i++) {
      if (bonded_ia_params[i].num != 2) {
        runtimeErrorMsg() << "The bonds for three particle binding need to be angle bonds.";
        return false;
      }
    }
  }
  
  // Create particle types


    if (collision_params.mode & COLLISION_MODE_VS) {
      if (n_nodes>1)
      {
        throw std::runtime_error("The make_particle_type_exists() logic in validate_collision_params() needs to be adapted once vs based methods support parallel simulations in the 1st place.");
      }
      if (collision_params.vs_particle_type<0){
        runtimeErrorMsg() << "Collision detection particle type for virtual sites needs to be >=0";
        return false;
      }
      make_particle_type_exist(collision_params.vs_particle_type);
    }

  
  
    if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)
    {
      if (n_nodes>1)
      {
        throw std::runtime_error("The make_particle_type_exists() logic in validate_collision_params() needs to be adapted once vs based methods support parallel simulations in the 1st place.");
      }
      if (collision_params.vs_particle_type<0){
        runtimeErrorMsg() << "Collision detection particle type for virtual sites needs to be >=0";
        return false;
      }
      make_particle_type_exist(collision_params.vs_particle_type);

      if (collision_params.part_type_to_be_glued<0){
        runtimeErrorMsg() << "Collision detection particle type to be glued needs to be >=0";
        return false;
      }
      make_particle_type_exist(collision_params.part_type_to_be_glued);
      
      if (collision_params.part_type_to_attach_vs_to<0){
        runtimeErrorMsg() << "Collision detection particle type to attach the virtual site to  needs to be >=0";
        return false;
      }
      make_particle_type_exist(collision_params.part_type_to_attach_vs_to);
      
      if (collision_params.part_type_after_glueing<0){
        runtimeErrorMsg() << "Collision detection particle type after glueing needs to be >=0";
        return false;
      }
      make_particle_type_exist(collision_params.part_type_after_glueing);
    }
  
  recalc_forces = 1;

  return true;
}

//* Allocate memory for the collision queue /
void prepare_collision_queue()
{
 TRACE(printf("%d: Prepare_collision_queue()\n",this_node));
  number_of_collisions=0;
}


inline bool bond_exists(const Particle* const p, const Particle* const partner, int bond_type)
{
  // First check the bonds of p1
  if (p->bl.e) {
    int i = 0;
    while(i < p->bl.n) {
      int size = bonded_ia_params[p->bl.e[i]].num;
      
      if (p->bl.e[i] == bond_type &&
          p->bl.e[i + 1] == partner->p.identity) {
        // There's a bond, already. Nothing to do for these particles
        return true;
      }
      i += size + 1;
    }
  }
  return false;
}


inline void queue_collision(int part1,int part2) {

    //Get memory for the new entry in the collision queue
    number_of_collisions++;
    if (number_of_collisions==1)
      collision_queue = (collision_struct *) malloc(number_of_collisions*sizeof(collision_struct));
    else
      collision_queue = (collision_struct *) realloc (collision_queue,number_of_collisions*sizeof(collision_struct));
    // Save the collision      
    collision_queue[number_of_collisions-1].pp1 = part1;
    collision_queue[number_of_collisions-1].pp2 = part2;
    
    TRACE(printf("%d: Added to queue: Particles %d and %d at %lf %lf %lf\n",this_node,part1,part2,point_of_collision[0],point_of_collision[1],point_of_collision[2]));
}

inline bool glue_to_surface_criterion(const Particle* const p1, const Particle* const p2) {
    return  (
       ((p1->p.type==collision_params.part_type_to_be_glued)
       && (p2->p.type ==collision_params.part_type_to_attach_vs_to))
      ||
       ((p2->p.type==collision_params.part_type_to_be_glued)
       && (p1->p.type ==collision_params.part_type_to_attach_vs_to)));
}

// Detect a collision between the given particles.
// Add it to the queue in case virtual sites should be added at the point of collision
void detect_collision(const Particle* const p1, const Particle* const p2, const double& dist_betw_part)
{

  if (dist_betw_part > collision_params.distance)
    return;

  //TRACE(printf("%d: particles %d and %d within bonding distance %lf\n", this_node, p1->p.identity, p2->p.identity, dist_betw_part));
  // If we are in the glue to surface mode, check that the particles
  // are of the right type
  if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)
    if (!glue_to_surface_criterion(p1,p2))
       return;
   

      
#ifdef VIRTUAL_SITES_RELATIVE
  // Ignore virtual particles
  if ((p1->p.isVirtual) || (p2->p.isVirtual))
    return;
#endif


  // Check, if there's already a bond between the particles
  if (bond_exists(p1,p2, collision_params.bond_centers))
    return;
  
  if (bond_exists(p2,p1, collision_params.bond_centers))
    return;


  TRACE(printf("%d: no previous bond, binding\n", this_node));

  /* If we're still here, there is no previous bond between the particles,
     we have a new collision */

  

  // do not create bond between ghost particles
  if (p1->l.ghost && p2->l.ghost) {
     TRACE(printf("Both particles %d and %d are ghost particles", p1->p.identity, p2->p.identity));
     return;
  }
  queue_collision(p1->p.identity,p2->p.identity);

  
}

/** @brief Calculate position of vs for GLUE_TO_SURFACE mode 
*    Reutnrs id of particle to bind vs to */
int glue_to_surface_calc_vs_pos(const Particle* const p1, const Particle* const p2, double pos[3]) 
{
    int bind_vs_to_pid;
    double vec21[3];
    double c;
    get_mi_vector(vec21,p1->r.p,p2->r.p);
    const double dist_betw_part=sqrt(sqrlen(vec21));
    
    // Find out, which is the particle to be glued.
    if ((p1->p.type==collision_params.part_type_to_be_glued)
       && (p2->p.type ==collision_params.part_type_to_attach_vs_to))
    { 
	       c = collision_params.dist_glued_part_to_vs/dist_betw_part;
         bind_vs_to_pid=p2->p.identity;
    }
    else if ((p2->p.type==collision_params.part_type_to_be_glued)
          && (p1->p.type ==collision_params.part_type_to_attach_vs_to))
    { 
	       c = -collision_params.dist_glued_part_to_vs/dist_betw_part;
         bind_vs_to_pid=p1->p.identity;
    }
    else
       {
         throw std::runtime_error("This should never be thrown. Bug.");
    }
    for (int i=0;i<3;i++) {
       pos[i] = p1->r.p[i] - vec21[i] * c;
    }
   return bind_vs_to_pid;
}

void bind_at_point_of_collision_calc_vs_pos(const Particle* const p1, const Particle* const p2, double pos1[3],double pos2[3]) {
    double vec21[3];
    get_mi_vector(vec21,p1->r.p,p2->r.p);
    for (int i=0;i<3;i++) {
       pos1[i] = p1->r.p[i] - vec21[i] * collision_params.vs_placement;
       pos2[i] = p1->r.p[i] - vec21[i] * (1.-collision_params.vs_placement);
    }
}





// Considers three particles for three_particle_binding and performs
// the binding if the criteria are met //
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

  //TRACE(printf("%d: proceeding to install three particle bond %d %d %d\n", this_node, p1->p.identity, p->p.identity, p2->p.identity));

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

// If activated, throws an exception for each collision which can be
// parsed by the script interface
void handle_exception_throwing_for_single_collision(int i)
{
    if (collision_params.exception_on_collision) {

      int id1, id2;
      if (collision_queue[i].pp1 > collision_queue[i].pp2) {
	id1 = collision_queue[i].pp2;
	id2 = collision_queue[i].pp1;
      }
      else {
	id1 = collision_queue[i].pp1;
	id2 = collision_queue[i].pp2;
      }
      std::ostringstream msg;
      msg << "collision between particles " << id1 << " and " <<id2;
      runtimeError(msg);
    }
}

#ifdef VIRTUAL_SITES_RELATIVE
void place_vs_and_relate_to_particle(double* pos, int relate_to)
{
 
	  place_particle(max_seen_particle+1,pos);
	  vs_relate_to(max_seen_particle,relate_to);
	  
	  (local_particles[max_seen_particle])->p.isVirtual=1;
	  #ifdef ROTATION_PER_PARTICLE
	    (local_particles[relate_to])->p.rotation=14;
	  #endif
	  (local_particles[max_seen_particle])->p.type=collision_params.vs_particle_type;
}


void bind_at_poc_create_bond_between_vs(int i)
{
   int bondG[3];

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

void glue_to_surface_bind_vs_to_pp1(int i)
{
	 int bondG[3];
         // Create bond between the virtual particles
         bondG[0] = collision_params.bond_vs;
         bondG[1] = max_seen_particle;
         local_change_bond(collision_queue[i].pp1, bondG, 0);
	 local_particles[collision_queue[i].pp1]->p.type=collision_params.part_type_after_glueing;
}

#endif

void gather_collision_queue(int* counts)
{
    int displacements[n_nodes];                   // offsets into collisions
  
    // Initialize number of collisions gathered from all processors
    for (int a=0;a<n_nodes;a++)
        counts[a]=0;
    
    // Total number of collisions
    MPI_Allreduce(&number_of_collisions, &total_collisions, 1, MPI_INT, MPI_SUM, comm_cart);
    
    if (total_collisions==0)
      return;

    // Gather number of collisions
    MPI_Allgather(&number_of_collisions, 1, MPI_INT, counts, 1, MPI_INT, comm_cart);

    // initialize displacement information for all nodes
    displacements[0]=0;
  
    // Find where to place collision information for each processor
    int byte_counts[n_nodes];
    for (int k=1; k<n_nodes; k++)
        displacements[k]=displacements[k-1]+(counts[k-1])*sizeof(collision_struct);
    
    for (int k=0; k<n_nodes; k++)
       byte_counts[k]=counts[k]*sizeof(collision_struct);
    
    TRACE(printf("counts [%d] = %d and number of collisions = %d and diplacements = %d and total collisions = %d\n", this_node, counts[this_node], number_of_collisions, displacements[this_node], total_collisions));
    
    // Allocate mem for the new collision info
    gathered_queue = (collision_struct *) malloc(total_collisions * sizeof(collision_struct));

    // Gather collision informtion from all nodes and send it to all nodes
    MPI_Allgatherv(collision_queue, byte_counts[this_node], MPI_BYTE, gathered_queue, byte_counts, displacements, MPI_BYTE, comm_cart);

    return;
}


// this looks in all local particles for a particle close to those in a 
// 2-particle collision. If it finds them, it performs three particle binding
void three_particle_binding_full_search()
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
}


// Goes through the collision queue and for each pair in it
// looks for a third particle by using the domain decomposition
// cell system. If found, it performs three particle binding
void three_particle_binding_domain_decomposition()
{
  // We have domain decomposition
    
  // Indices of the cells in which the colliding particles reside
  int cellIdx[2][3];
    
  // Iterate over collision queue

  for (int id=0;id<total_collisions;id++) {

      // Get first cell Idx
      if ((local_particles[gathered_queue[id].pp1] != NULL) && (local_particles[gathered_queue[id].pp2] != NULL)) {

        Particle* p1=local_particles[gathered_queue[id].pp1];
        Particle* p2=local_particles[gathered_queue[id].pp2];
        dd_position_to_cell_indices(p1->r.p,cellIdx[0]);
        dd_position_to_cell_indices(p2->r.p,cellIdx[1]);

        // Iterate over the cells + their neighbors
        // if p1 and p2 are in the same cell, we don't need to consider it 2x
        int lim=1;

        if ((cellIdx[0][0]==cellIdx[1][0]) && (cellIdx[0][1]==cellIdx[1][1]) && (cellIdx[0][2]==cellIdx[1][2]))
          lim=0; // Only consider the 1st cell

        for (int j=0;j<=lim;j++) {

            // Iterate the cell with indices cellIdx[j][] and all its neighbors.
            // code taken from dd_init_cell_interactions()
            for(int p=cellIdx[j][0]-1; p<=cellIdx[j][0]+1; p++)	
               for(int q=cellIdx[j][1]-1; q<=cellIdx[j][1]+1; q++)
	                for(int r=cellIdx[j][2]-1; r<=cellIdx[j][2]+1; r++) {   
	                   int ind2 = get_linear_index(p,q,r,dd.ghost_cell_grid);
	                   Cell* cell=&cells[ind2];
 
	                   // Iterate over particles in this cell
                     for(int a=0; a<cell->n; a++) {
                        Particle* P=&cell->part[a];
                        // for all p:
  	                      // Check, whether p is equal to one of the particles in the
  	                      // collision. If so, skip
  	                      if ((P->p.identity ==p1->p.identity) || (P->p.identity == p2->p.identity)) {
                          //TRACE(printf("same particle\n"));
  		                continue;
  	                      }

                        // The following checks, 
                        // if the particle p is closer that the cutoff from p1 and/or p2.
                        // If yes, three particle bonds are created on all particles
  	                      // which have two other particles within the cutoff distance,
                        // unless such a bond already exists

  	                      // We need all cyclical permutations, here 
                        // (bond is placed on 1st particle, order of bond partners
  	                      // does not matter, so we don't need non-cyclic permutations):

                        if (P->l.ghost) {
                          //TRACE(printf("%d: center particle is ghost: %d\n", this_node, P->p.identity));
                          continue;
                        }
                        //TRACE(printf("%d: LOOP: %d Handling collision of particles FIRST CONFIGURATION %d %d %d\n", this_node, id, p1->p.identity, P->p.identity, p2->p.identity));
                        coldet_do_three_particle_bond(P,p1,p2);

                        if (p1->l.ghost) {
                          //TRACE(printf("%d: center particle is ghost: %d\n", this_node, p1->p.identity));
                          continue;
                        }

                        coldet_do_three_particle_bond(p1,P,p2);

                        if (p2->l.ghost) {
                          //TRACE(printf("%d: center particle is ghost: %d\n", this_node, p2->p.identity));
                          continue;
                        }

  	                      coldet_do_three_particle_bond(p2,P,p1);

                     } // loop over particles in this cell

	                } // Loop over cell

        } // Loop over particles if they are in different cells
    
      } // If local particles exist

  } // Loop over total collisions
}


// Handle the collisions stored in the queue
void handle_collisions ()
{

  TRACE(printf("%d: handle_collisions: number of collisions in queue %d\n",this_node,number_of_collisions));  

  if (collision_params.exception_on_collision) {
    for (int i=0;i<number_of_collisions;i++) {
      handle_exception_throwing_for_single_collision(i);
    } 
  }
    
    
  if (bind_centers()) 
  {
    for (int i=0;i<number_of_collisions;i++) {
      // put the bond to the physical particle; at least one partner always is
      int primary =collision_queue[i].pp1;
      int secondary = collision_queue[i].pp2;
      if (local_particles[collision_queue[i].pp1]->l.ghost) {
        primary = collision_queue[i].pp2;
        secondary = collision_queue[i].pp1;
        TRACE(printf("%d: particle-%d is ghost", this_node, collision_queue[i].pp1));
      }
      int bondG[2];
      bondG[0]=collision_params.bond_centers;
      bondG[1]=secondary;
      local_change_bond(primary, bondG, 0);
      TRACE(printf("%d: Adding bond %d->%d\n",this_node, primary,secondary));
    }
  }

#ifdef VIRTUAL_SITES_RELATIVE
  // If one of the collision modes is active which places virtual sites, we go over the queue to handle them
  if ((collision_params.mode & COLLISION_MODE_VS) || (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)) {
    for (int i=0;i<number_of_collisions;i++) {
	// Create virtual site(s) 
  const int primary=collision_queue[i].pp1;
  const int secondary=collision_queue[i].pp2;
  const Particle* const p1=local_particles[primary];
  const Particle* const p2=local_particles[secondary];
	
	// If we are in the two vs mode
	// Virtual site related to first particle in the collision
	if (collision_params.mode & COLLISION_MODE_VS)
	{
   double pos1[3],pos2[3];

   bind_at_point_of_collision_calc_vs_pos(p1,p2,pos1,pos2);
	 place_vs_and_relate_to_particle(pos1,primary);
	 place_vs_and_relate_to_particle(pos2,secondary);
   bind_at_poc_create_bond_between_vs(i);
	}
	if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF) {
      double pos[3];
      const int pid=glue_to_surface_calc_vs_pos(p1,p2,pos);
      place_vs_and_relate_to_particle(pos,pid);
      glue_to_surface_bind_vs_to_pp1(i);
    }
      } // Loop over all collisions in the queue
    } // are we in one of the vs_based methods
#endif //defined VIRTUAL_SITES_RELATIVE
  

  // three-particle-binding part


  if (collision_params.mode & (COLLISION_MODE_BIND_THREE_PARTICLES)) {
  int counts[n_nodes];
  gather_collision_queue(counts);

    if (counts[this_node]>0) {

      // If we don't have domain decomposition, we need to do a full sweep over all
      // particles in the system. (slow)
      if (cell_structure.type!=CELL_STRUCTURE_DOMDEC) {
        three_particle_binding_full_search();
    } // if cell structure != domain decomposition
    else
    {
      three_particle_binding_domain_decomposition();
    } // If we have doamin decomposition

   } // if number of collisions of this node > 0
       
       if (total_collisions>0)
         free(gathered_queue);
       total_collisions = 0;
 } // if TPB

  // If a collision method is active which places particles, resorting might be needed
  TRACE(printf("%d: Resort particles is %d\n",this_node,resort_particles));
  if (collision_params.mode & (COLLISION_MODE_VS | COLLISION_MODE_GLUE_TO_SURF))
  {
    // NOTE!! this has to be changed to total_collisions, once parallelization
    // is implemented

    if (number_of_collisions >0)
    {
      announce_resort_particles();
    }
  }
  
  // Reset the collision queue
  if (number_of_collisions>0)
    free(collision_queue);
  
  number_of_collisions = 0;


}

#endif
