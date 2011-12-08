#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "collision.h"
#include "virtual_sites.h"
#include "integrate.h"
#include "cells.h"
#include "communication.h" 
#include "parser.h" 


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
extern int collision_detection_mode=0;

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

// Set params
collision_detection_mode=mode;
collision_detection_bond_centers=bond_centers;
collision_detection_bond_vs=bond_vs;
collision_distance=d;
collision_vs_particle_type=t;
make_particle_type_exist(t);
recalc_forces=1;
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
  int bondG[2], i;

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
  
      // Create bond between the virtual particles
      bondG[0] =collision_detection_bond_vs;
      bondG[1] =max_seen_particle-1;
      local_change_bond(max_seen_particle, bondG, 0);
    }

    // Is this needed
    //on_particle_change();
  }

}
#endif

  // Reset the collision queue	 
  number_of_collisions = 0;
  free(collision_queue);

}


int tclcommand_on_collision(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
 // If no argumens are given, print status
 if (argc==1)
 {
  char s[200];
  if (collision_detection_mode==0)
  {
   sprintf(s,"off");
  }
  
  if (collision_detection_mode==1)
  {
   sprintf(s,"bind_centers %f %d",collision_distance,collision_detection_bond_centers);
  }
  
  if (collision_detection_mode==2)
  {
   sprintf(s,"bind_at_point_of_collision %f %d %d %d",collision_distance,collision_detection_bond_centers, collision_detection_bond_vs,collision_vs_particle_type);
  }
  Tcl_AppendResult(interp, s, (char*) NULL);
  return TCL_OK;
 }

 // Otherwise, we set parametes

 int res;

 if (ARG1_IS_S("off"))
 {
  collision_detection_set_params(0,0.,0,0,0);
 }
 else if (ARG1_IS_S("bind_centers"))
 {
  if (argc!=4)
  {
   Tcl_AppendResult(interp, "Need a ditance and a bond type as args.", (char*) NULL);
   return TCL_ERROR;
  }
  double d;
  if (!ARG_IS_D(2,d))
  {
   Tcl_AppendResult(interp, "Need a ditance as 1st arg.", (char*) NULL);
   return TCL_ERROR;
  }
  int bond1;
  if (!ARG_IS_I(3,bond1))
  {
   Tcl_AppendResult(interp, "Need a bond type as 2nd argument.", (char*) NULL);
   return TCL_ERROR;
  }
  res=collision_detection_set_params(1,d,bond1,0,0);
  if (res==2)
  {
   Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==3)
  {
   Tcl_AppendResult(interp, "A bond of the specified type does not exist.", (char*) NULL);
   return TCL_ERROR;
  }
 }
 else if (ARG1_IS_S("bind_at_point_of_collision"))
 {
  if (argc!=6)
  {
   Tcl_AppendResult(interp, "Need a ditance, two bond types, and a particle type as args.", (char*) NULL);
   return TCL_ERROR;
  }
  double d;
  if (!ARG_IS_D(2,d))
  {
   Tcl_AppendResult(interp, "Need a ditance as 1st arg.", (char*) NULL);
   return TCL_ERROR;
  }
  int bond1,bond2,t;
  if ((!ARG_IS_I(3,bond1)) || (!ARG_IS_I(4,bond2)) || (!ARG_IS_I(5,t)))
  {
   Tcl_AppendResult(interp, "Need two bond types as 2nd and 3rd and a particle type as 4th argument.", (char*) NULL);
   return TCL_ERROR;
  }
  res=collision_detection_set_params(2,d,bond1,bond2,t);
  if (res==1)
  {
   Tcl_AppendResult(interp, "This mode requires the VIRTUAL_SITES_RELATIVE feature to be comiled in.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==2)
  {
   Tcl_AppendResult(interp, "Collision detection only works on a single cpu.", (char*) NULL);
   return TCL_ERROR;
  }
  if (res==3)
  {
   Tcl_AppendResult(interp, "A bond of the specified type does not exist.", (char*) NULL);
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
