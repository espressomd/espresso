/**************************************************/
/*******************  FORCES.C  *******************/
/**************************************************/

#include "forces.h"
#include "debug.h"

void force_init()
{

  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d interaction types\n",
		      this_node,n_particle_types));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_interaction_types));
}

void force_calc()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_calc:\n",this_node));
}

void force_exit()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_exit:\n",this_node));
}
