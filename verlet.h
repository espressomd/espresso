/**************************************************/
/*******************  VERLET.H  *******************/
/**************************************************/
#ifndef VERLET_H
#define VERLET_H
#include "global.h"

/** initialize verlet list structure. */
void verlet_init();

/** fill the verlet table. */
void build_verlet_list();

#endif
