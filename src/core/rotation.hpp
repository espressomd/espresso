/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file rotation.hpp
    This file contains all subroutines required to process rotational motion.

*/

#include "utils.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"
#include "gb.hpp"

/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/
 
/** Propagate angular velocities and update quaternions on a particle */
void propagate_omega_quat_particle(Particle* p); 

/** Convert torques to the body-fixed frame and propogate
    angular velocities */
void convert_torques_propagate_omega();

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques();

/** convert angular velocities and torques from the 
    body-fixed frames to space-fixed coordinates */
void convert_omega_body_to_space(Particle *p, double *omega);
void convert_torques_body_to_space(Particle *p, double torque[3]);

/** Here we use quaternions to calculate the rotation matrix which
    will be used then to transform torques from the laboratory to
    the body-fixed frames */  
void define_rotation_matrix(Particle *p, double A[9]);

inline void convert_quat_to_quatu(double quat[4], double quatu[3])
{
  /* director */
  quatu[0] = 2*(quat[1]*quat[3] + quat[0]*quat[2]);
  quatu[1] = 2*(quat[2]*quat[3] - quat[0]*quat[1]);
  quatu[2] =   (quat[0]*quat[0] - quat[1]*quat[1] -
		quat[2]*quat[2] + quat[3]*quat[3]); 
}

/** Multiply two quaternions */ 
void multiply_quaternions(double a[4], double b[4], double result[4]);

/** Convert director to quaternions */
int convert_quatu_to_quat(double d[3], double quat[4]);


#ifdef DIPOLES

/** convert a dipole moment to quaternions and dipolar strength  */
inline int convert_dip_to_quat(double dip[3], double quat[4], double *dipm)
{
  double dm;
  // Calculate magnitude of dipole moment
  dm = sqrt(dip[0]*dip[0] + dip[1]*dip[1] + dip[2]*dip[2]);
  *dipm = dm;
  convert_quatu_to_quat(dip,quat);

  return 0;
}

/** convert quaternion director to the dipole moment */
inline void convert_quatu_to_dip(double quatu[3], double dipm, double dip[3])
{
  /* dipole moment */
  dip[0] = quatu[0]*dipm;
  dip[1] = quatu[1]*dipm;
  dip[2] = quatu[2]*dipm;
}

#endif

#endif
