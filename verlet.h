/**************************************************/
/*******************  VERLET.H  *******************/
/**************************************************/
#ifndef VERLET_H
#define VERLET_H

#include <stdio.h>
#include <stdlib.h>
#include "global.h"



/** initialize verlet list structure. */
void verlet_init();

/** fill the verlet table. */
void build_verlet_list();

/** exit verlet list structure. */
void verlet_exit();

#endif
