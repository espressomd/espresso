// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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
