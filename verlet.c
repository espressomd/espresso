/**************************************************/
/*******************  VERLET.C  *******************/
/**************************************************/

#include "verlet.h"
#include "debug.h"

int rebuild_verletlist = 1;

void verlet_init()
{
  VERLET_TRACE(fprintf(stderr,"%d: verlet_init:\n",this_node));
}

void build_verlet_list()
{
  VERLET_TRACE(fprintf(stderr,"%d: build_verlet_list:\n",this_node));

  rebuild_verletlist = 0;
}

void verlet_exit()
{
  VERLET_TRACE(fprintf(stderr,"%d: verlet_exit:\n",this_node));
}
