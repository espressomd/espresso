#ifndef ROTATION_H
#define ROTATION_H
/** \file rotation.h
    This file ontains all subroutines required to process rotational motion. 

    <b>Responsible:</b>
    <a href="mailto:antypov@mpip-mainz.mpg.de">Dmytro</a>

*/

#include <tcl.h>
#include "config.h"
#include "particle_data.h"
#include "thermostat.h"
#include "gb.h"

/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/
 
/** Propagate angular velocities and update quaternions */
void propagate_omega_quat(); 

/** Convert torques to the body-fixed frame and propogate
    angular velocity */
void convert_torqes_propagate_omega();

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques();

#endif
