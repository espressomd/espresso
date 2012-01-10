#ifndef COLLISION_H
#define COLLISION_H
//#ifdef COLLISION_DETECTION


#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "virtual_sites.h"
#include "integrate.h"


// Data teyp holding the info about a single collision

typedef struct {
int pp1; // 1st particle id
int pp2; // 2nd particle id
double point_of_collision[3]; 
} collision_struct;


// Collision detection mod
// 0=off
// 1=bind centers
// 2=bind at point of collision
// Never write to this variable. Use collision_detection_set_params()
extern int collision_detection_mode;





// Detect a collision between two particles. In case of collision,
// a bond between the particles is added as marker and the collision is
// recorded in the queue for later processing.
void detect_collision(Particle* p1, Particle* p2);

void prepare_collision_queue();


// Handle the collisions recorded in the queue
void handle_collisions();


int tclcommand_on_collision(ClientData data, Tcl_Interp *interp, int argc, char **argv); 
//#endif
#endif
