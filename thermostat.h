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

/************************************************
 * exported variables
 ************************************************/

extern double friction_gamma;
extern double temperature;

/************************************************
 * functions
 ************************************************/

/** overwrite the forces of the local particles with
    the friction term, i. e. \f$ F_i=-\gamma v_i\f$.
*/
void friction_thermo();

/** Callback for setting \ref temperature */
int temp_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref friction_gamma */
int gamma_callback(Tcl_Interp *interp, void *_data);

#endif
