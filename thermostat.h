#ifndef THERMOSTAT_H
#define THERMOSTAT_H
/** \file thermostat.h
    Contains all thermostats. 

    <b>Responsible:</b>
    <a href="mailto:muehlbac@mpip-mainz.mpg.de">Frank</a>


    Currently there is only the
    \ref friction_thermo.
*/
#include <tcl.h>
#include "particle_data.h"

/************************************************
 * exported variables
 ************************************************/

extern double friction_gamma;
extern double temperature;

/************************************************
 * functions
 ************************************************/

/** initialize constants of the thermostat on
    start of integration */
void thermo_init();

/** overwrite the forces of a particle with
    the friction term, i.e. \f$ F_i= -\gamma v_i + \xi_i\f$.
*/
void friction_thermo(Particle *p);

/** set the particle torques to the friction term, i.e. \f$\tau_i=-\gamma w_i + \xi_i\f$.

The same friction coefficient \f$\gamma\f$ is used as that for translation.
*/
void friction_thermo_rotation(Particle *p);

/** Callback for setting \ref #temperature */
int temp_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref friction_gamma */
int gamma_callback(Tcl_Interp *interp, void *_data);

#endif
