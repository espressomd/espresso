/**************************************************/
/*******************  FORCES.C  *******************/
/**************************************************/

#include "forces.h"

#define DEBUG

void force_init()
{

#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: force_init:\n",this_node);
  if(this_node<2) fprintf(stderr,"%d: found %d interaction types\n",
			  this_node,n_particle_types);
  if(this_node<2) fprintf(stderr,"%d: found %d particles types\n",
			  this_node,n_interaction_types);
#endif

}

void force_calc()
{
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: force_calc:\n",this_node); 
#endif
}

void force_exit()
{
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: force_exit:\n",this_node); 
#endif
}
