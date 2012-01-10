#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "collision.h"
#include "virtual_sites.h"
#include "integrate.h"
#include "cells.h"


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

// Bonding mode: 
// 0 = off
// 1 = only create bond between centers of colliiding particles
// 2 = also create two virtual sites at the point o collision and bind the 
// together. This prevents the particles from sliding against each other. 
// (Requires VIRTUAL_SITES_RELATIVE)
int collision_detection_mode=0;

int seet_collision_detection_params(int mode, int bond_centers, int, bond_vs)
{
// If we don't have virtual sites relatie feature, mode 2 isn't allowed.
#ifndef VIRTUAL_SITES_RELATIVE
if (mode==2)
 return 1;
}
#endif

// Collision detection only works on a single cpu
if (mode!=0 && n_nodes !=1)
 return 2;

// Set params
collision_detection_mode=mode;
collision_detection_bond_centers=bond_centers;
collision_detection_bond_vs=bond_vs;
}

// Detect a collision between the given particles.
// In case of a collision, a bond is added between them as marker
// and the collision is recorded in the queue
void detect_collision(Particle* p1, Particle* p2)
{
  //printf("in collsiion_detction"); 
  double dist_betw_part, vec21[3], collisioncriter=1.15;
  int part1, part2, the_bond_type_added_on_collision=0, size;

  // Obtain distance between particles
  dist_betw_part = distance2vec(p1->r.p, p2->r.p, vec21);
  if (dist_betw_part > collisioncriter)
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
  int found = 0;
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

  // Create the marker bond between the particles
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
    //printf("connection point is %f\n",new_position[i]);
  }
       
  number_of_collisions = number_of_collisions+1;
  //	printf("number of collisions each time are %d\n",number_of_collisions);       
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
// Does the user wnat bonds between virtual sites placed at the point of collision
if (collision_detection_mode==2)
{

  //	printf("number of collisions in handle collision are %d\n",number_of_collisions);  
  int delete =0, bondG[2], i;

  if (number_of_collisions > 0) {
    // Go through the queue
    for (i=0;i<number_of_collisions;i++) {
      //  printf("Handling collision of particles %d %d\n", collision_queue[i].pp1, collision_queue[i].pp2);
      //  fflush(stdout);
   
      int j;

      // The following lines will remove the relative velocity from
      // colliding particles
      //   double v[3];
      //   for (j=0;j<3;j++)
      //   {
      //    v[j] =0.5 *((local_particles[collision_queue[i].pp1])->m.v[j] +(local_particles[collision_queue[i].pp2])->m.v[j]);
      //    (local_particles[collision_queue[i].pp1])->m.v[j] =v[j];
      //    (local_particles[collision_queue[i].pp2])->m.v[j] =v[j];
      //   }

      // Create virtual sites and bind them gotether
  
      // Virtual site related to first particle in the collision
      place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
      vs_relate_to(max_seen_particle,collision_queue[i].pp1);
      (local_particles[max_seen_particle])->p.isVirtual=1;
  
      //printf("virtual1 was created: %d rel_to %d\n", max_seen_particle, collision_queue[i].pp1); 
      // Virtual particle related to 2nd particle of the collision
      place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
      vs_relate_to(max_seen_particle,collision_queue[i].pp2);
      (local_particles[max_seen_particle])->p.isVirtual=1;
      //printf("virtual2 was created: %d rel_to %d\n", max_seen_particle, collision_queue[i].pp2); 
  
      // Create bond between the virtual particles
      bondG[0] =collision_detection_bond_marker_virtual;
      bondG[1] =max_seen_particle-1;
      local_change_bond(max_seen_particle, bondG, 0);
      //printf("virtual bond was created\n");   
    }

    on_particle_change();
  }

  // Reset the collision queue	 
}
#endif

  number_of_collisions = 0;
  free(collision_queue);

}


int tclcommand_on_collision(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
 // If no argumens are given, print status
 if (argc==1)
 {
  char[200] s;
  if (mode==0)
  {
   sprintf(s,"off");
  }
  
  if (mode==1)
  {
   sprintf(s,"bind_centers %d",collision_detection_bond_centers);
  }
  
  if (mode==2)
  {
   sprintf(s,"bind_at_point_of_collision %d %d",collision_detection_bond_centers, collision_detection_bond_vs);
  }
 }

 // Otherwise, we set parametes
 char[80] modestr;
 ARG1_IS_S(modestr);

 int res;

 if (modestr=="off")
 {
  collision_detection_set_params(mode,0,0);
 }
 else if (modestr=="bind_centers")
 {
  int bond1;
  if (!ARG_IS_I(2,bond1)
  {
   Tcl_AppendResult(interp, "Need a bond type as argument.", (char*) NULL);
   return TCL_ERROR;
  }
  res=collision_detection_set_params(mode,bond1,0);
  if (res==2)
  {
   Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
   return TCL_ERROR;
  }
 }
 if (modestr=="bind_at_point_of_collision")
 {
  int bond1,bond2;
  if ((!ARG_IS_I(2,bond1) || (!ARG_IS_I(3,bond2)))
  {
   Tcl_AppendResult(interp, "Need two bond types as argument.", (char*) NULL);
   return TCL_ERROR;
  }
  res=collision_detection_set_params(mode,bond1,bond1);
  if (res==1)
  {
   Tcl_AppendResult(interp, "This mode requires the VIRTUAL_SITES_RELATIVE feature to be complied in.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==2)
  {
   Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
   return TCL_ERROR;
  }
 }
 else
 {
  Tcl_AppendResult(interp,"Unknown mode.",(char*)NULL);
  return TCL_ERROR;
 }
 return TCL_OK;
}




#endif
