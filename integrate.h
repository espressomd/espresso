#ifndef INTEGRATE_H
#define INTEGRATE_H
/** \file integrate.h    Molecular dynamics integrator.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  For more information see \ref integrate.c "integrate.c".
*/   
#include <tcl.h>

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Time step for the integration. */
extern double time_step;
/** Old time step needed for rescaling of forces. */
extern double old_time_step;
/** Physical start time of the simulation. */
extern double start_time;
/** Actual simulation time (only on MASTER NODE). */
extern double sim_time;
/** Maximal interaction cutoff. */
extern double max_cut;
/** Verlet list skin. */
extern double skin;
/** Maximal interaction range (max_cut + skin). */
extern double max_range;
/** square of \ref max_range. */
extern double max_range2;
/** If non-zero, some particles have moved since last
    integration. */
extern int    particle_changed;
/** If non-zero, some non-bonded interactions have changed since last
    integration. */
extern int    interactions_changed;
/** If non-zero, the system topology has changed since last
    integration. */
extern int    topology_changed;
/** If non-zero, some other parameter (e.g. time_step, skin) has changed
    since last integration. */
extern int    parameter_changed;
/** Average number of integration steps the verlet list has been re
    used. */
extern double verlet_reuse;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure for integrator steering.
    USAGE: integrate <steps> \\   
    see also \ref tcl_integrate
*/
int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** Calculate maximal interaction range. 
    Uses \ref calc_maximal_cutoff.
    \ref max_range  = \ref max_cut + \ref #skin;
 */
void integrate_vv_recalc_maxrange();

/** integrate with velocity verlet integrator.
    \param n_steps number of steps to integrate.
 */
void integrate_vv(int n_steps);

/** function that rescales all velocities on one node according to a
    new time step. */
void rescale_velocities(); 

/** Callback for setmd skin.
    \return TCL status.
*/
int skin_callback(Tcl_Interp *interp, void *_data);

/** Callback for integration time_step (0.0 <= time_step).
    \return TCL status.
*/
int time_step_callback(Tcl_Interp *interp, void *_data);

/** Callback for start_time of the integration.
    If no value is set the integration starts at start_time = 0.0.
    \return TCL status.
*/
int start_time_callback(Tcl_Interp *interp, void *_data);

/*@}*/

#endif
