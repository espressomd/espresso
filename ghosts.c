/**************************************************/
/*******************  GHOSTS.C  *******************/
/**************************************************/

#include "ghosts.h"
#include "debug.h"

/** flag for double sided ghost frame. */
int double_sided = 0;

void ghost_init()
{
  GHOST_TRACE(fprintf(stderr,"%d: ghost_init:\n",this_node));

  // if(n_cells==0) {
  //  fprintf(stderr,"%d Cell structure not initialized! STOP!",this_node);
  //  exit(1);
  // }
}

void exchange_part()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_part:\n",this_node));
}

void exchange_ghost()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost:\n",this_node));
}

void exchange_ghost_pos()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost_pos:\n",this_node));
}

void exchange_ghost_forces()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost_forces:\n",this_node));
}

void ghost_exit()
{
  GHOST_TRACE(fprintf(stderr,"%d: ghost_exit:\n",this_node));
}
