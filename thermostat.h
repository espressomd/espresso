#ifndef THERMOSTAT_H
#define THERMOSTAT_H
/** \file thermostat.h
    Contains all thermostats. Currently there is only the
    \ref friction_thermo.
*/

/************************************************
 * exported variables
 ************************************************/

extern double friction_gamma;

/************************************************
 * functions
 ************************************************/

/** overwrite the forces of the local particles with
    the friction term, i. e. \f$ F_i=-\gamma v_i\f$.
*/
void friction_thermo();

#endif
