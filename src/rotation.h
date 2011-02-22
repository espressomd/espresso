/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef ROTATION_H
#define ROTATION_H
/** \file rotation.h
    This file contains all subroutines required to process rotational motion.

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
void convert_torques_propagate_omega();

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques();

/** convert torques from the body-fixed frames to space-fixed coordinates */
void convert_torques_body_to_space(Particle *p, double torque[3]);

/** convert quaternions to the director */
void convert_quat_to_quatu(double quat[4], double quatu[3]);

/** convert dipole moment of one particle to the quaternions. Returns 1 if
    everything is ok, or 0 if the dipole vector was too small. */
int convert_dip_to_quat(double dip[3], double quat[4], double *dipm);

/** convert quaternion director to the dipole moment */
void convert_quatu_to_dip(double quatu[3], double dipm, double dip[3]);

/** Multiply two quaternions */ 
void multiply_quaternions(double a[4], double b[4], double result[4]);

/** Convert director to quaternions */
int convert_quatu_to_quat(double d[3], double quat[4]);

void convert_omega_body_to_space(Particle *p, double *omega);
#endif
