/**************************************************/
/*******************  VERLET.C  *******************/
/**************************************************/

#include "verlet.h"

//#define DEBUG

int rebuild_verletlist = 1;

void verlet_init()
{
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: verlet_init:\n",this_node); 
#endif
}

void build_verlet_list()
{
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: build_verlet_list:\n",this_node); 
#endif
  rebuild_verletlist = 0;
}

void verlet_exit()
{
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: verlet_exit:\n",this_node); 
#endif
}
