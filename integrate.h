/*****************************************************/
/*******************  INTEGRATE.H  *******************/
/*****************************************************/
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include <math.h>
#include "communication.h"
#include "global.h"
#include "cells.h"
#include "verlet.h"
#include "forces.h"
#include "ghosts.h"

/** declare commands needed for steering the integrator. */
void tcl_integrator_init(Tcl_Interp *interp);

/** tcl procedure for integrator steering */
int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** initialize velocity verlet integrator. */
void integrate_vv_init();

/** integrate with velocity verlet integrator. */
void integrate_vv(int n_steps);

/** exit velocity verlet integrator. */
void integrate_vv_exit();

#endif
