// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef INTEGRATE_H
#define INTEGRATE_H
/** \file integrate.h    Molecular dynamics integrator.
 *
 *  For more information see \ref integrate.c "integrate.c".
*/   
#include <tcl.h>

#define INTEG_METHOD_NPT_ISO   0
#define INTEG_METHOD_NVT       1

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Switch determining which Integrator to use. */
extern int integ_switch;

/** incremented if a Verlet update is done, aka particle resorting. */
extern int n_verlet_updates;

/** Time step for the integration. */
extern double time_step;
/** Old time step needed for rescaling of forces. */
extern double old_time_step;
/** Actual simulation time (only on MASTER NODE). */
extern double sim_time;
/** Maximal interaction cutoff. */
extern double max_cut;
/** Verlet list skin. */
extern double skin;
/** Maximal interaction range (max_cut + skin). */
extern double max_range;
/** Square of \ref max_range. It's initial value is -1.0 which is
    used to determine wether max_range/max_range2 has been set
    properly by \ref integrate_vv_recalc_maxrange or not. */
extern double max_range2;
/** Maximal non bonded interaction range (max_cut_non_bonded + skin). */
extern double max_range_non_bonded;
/** Square of \ref max_range_non_bonded. */
extern double max_range_non_bonded2;


/** If non-zero, the particle data will be resorted before the next integration. */
extern int    resort_particles;
/** If non-zero, the forces will be recalculated before the next integration. */
extern int    recalc_forces;
/** Average number of integration steps the verlet list has been re
    used. */
extern double verlet_reuse;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** tcl procedure for integrator steering. For documentation,
    see \ref tcl_integrate
*/
int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);

/** Calculate maximal interaction range. 
    Uses \ref calc_maximal_cutoff.
    \ref max_range  = \ref max_cut + \ref #skin;
 */
void integrate_vv_recalc_maxrange();

/** Initialize the used thermodynamic Ensemble (NVT or NPT) */
void integrate_ensemble_init();

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

/** Callback for current time in the integration.
    If no value is set the integration starts at time = 0.0.
    \return TCL status.
*/
int time_callback(Tcl_Interp *interp, void *_data);

/** Implements the tcl-command 'invalidate_system' which forces a system re-init. 
    For more information, see \ref tcl_invalidate_system. */
int invalidate_system(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** local routine of \ref invalidate_system */
void local_invalidate_system();

/*@}*/

#endif
