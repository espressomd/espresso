/*****************************************************/
/*******************  INTEGRATE.H  *******************/
/*****************************************************/
#ifndef INTEGRATE_H
#define INTEGRATE_H
#include <tcl.h>

/****************************************
 * exported variables
 ****************************************/

/** Time step for the integration. */
extern double time_step;
/** Maximal interaction cutoff. */
extern double max_cut;
/** Verlet list skin. */
extern double skin;
/** Maximal interaction range (max_cut + skin). */
extern double max_range;
/** square of \ref max_range. */
extern double max_range2;
/** Flag for integrator. If non-zero, the forces have to be calculated
    before the first step. */
extern int calc_forces_first;

/**********************************************
 * functions
 **********************************************/

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
