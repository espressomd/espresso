/**************************************************/
/*******************  FORCES.C  *******************/
/**************************************************/
#ifndef FORCES_H
#define FORCES_H

#include "forces.h"

#define DEBUG

void force_init()
{
  int i,j;
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: force_init:\n",this_node);
  if(this_node<2) fprintf(stderr,"%d: found %d interaction types\n",
			  this_node,n_particle_types);
  if(this_node<2) fprintf(stderr,"%d: found %d particles types\n",
			  this_node,n_interaction_types);
#endif

  

}

#endif
