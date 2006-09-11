// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef ROTATION_H
#define ROTATION_H
/** \file rotation.h
    This file contains all subroutines required to process rotational motion.

    <b>Responsible:</b>
    <a href="mailto:antypov@mpip-mainz.mpg.de">Dmytro</a>

*/

#include <tcl.h>
#include "utils.h"
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
    angular velocities */
void convert_torqes_propagate_omega();

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques();

/** convert torques from the body-fixed frames to space-fixed coordinates */
void convert_torques_body_to_space(Particle *p, double torque[3]);

/** convert omega to the body-fixed frames */
void convert_initial_omega(Particle *p, double omega[3]);

/** convert omega from the body-fixed frames to space-fixed coordinates */
void convert_omega_body_to_space(Particle *p, double omega[3]);

/** convert dipole moment of one particle to the quaternions  */
int convert_dip_to_quat_one(double dip[3], double quat[4]);

/** convert dipole moments of particles to the quaternions  */
void convert_dip_to_quat_all();

/** convert quaternions to the dipole moment of the particle  */
void convert_quat_to_dip_one(double dip[3], double dipm, double quat[4]);

/** convert quaternions to the dipole moments of particles  */
void convert_quat_to_dip_all();


#endif
