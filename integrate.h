/*******************  INTEGRATE.H *******************/
#ifndef INTEGRATE_H
#define INTEGRATE_H
#include "global.h"
#include <tcl.h>

/** declare commands needed for steering the integrator. */
void tcl_integrator_init(Tcl_Interp *interp);

/** tcl procedure for integrator steering */
int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** initialize velocity verlet integrator. */
int integrate_vv_init();

/** integrate with velocity verlet integrator. */
int integrate_vv(int n_steps);

/** exit velocity verlet integrator. */
int integrate_vv_exit();

#endif
