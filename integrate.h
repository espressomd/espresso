/*****************************************************/
/*******************  INTEGRATE.H  *******************/
/*****************************************************/
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <math.h>
#include "communication.h"
#include "global.h"
#include "cells.h"
#include "verlet.h"
#include "ghosts.h"
#include "forces.h"

/** tcl procedure for integrator steering.
    USAGE: integrate <task>                       \\
           task can either be init, #steps, exit. \\
    EXAMPLE for an integration:                   \\
           integrate init                         \\
           integrate 100                          \\
           integrate exit                         
*/
int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** initialize velocity verlet integrator. */
void integrate_vv_init();

/** integrate with velocity verlet integrator. */
void integrate_vv(int n_steps);

/** exit velocity verlet integrator. */
void integrate_vv_exit();

#endif
