/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
/** \file ibm.hpp
 * Header file for ibm.cpp
 *
 * This is the header file for some of the Immersed boundary method functions
 */

#ifndef IBM_H
#define IBM_H

#ifdef IMMERSED_BOUNDARY

/** Calculates the coupling of virtual_sites_ibm to the LB fluid.
 * This function  is called from \ref force_calc. The force is retrieved
 * from the particle force and then it is distributed to the fluid
 * Note that this function changes the state of the fluid!
 */
void lb_ibm_coupling();

/** Fluid-Membrane Interaction distributes particle force to surrounding fluid nodes via linear interpolation.
 * 
 * Note: At the time this function is called, the forces saved in p are not yet scaled! Only called for virtual parts
 *
 * @param p          The coupled particle (Input).
 */
void couple_trace_to_fluid(Particle *p);

/** Fluid Membrane interaction, calculate interpolated velocity with updated (+f_ext), but not yet streamed modes
 * 
 * follows the implementation of lb_fluid_get_interpolated_velocity in lb.cpp but with forces on nodes 
 * taken into account by linear interpolation
 *
 * @param p          The coupled particle (Input).
 */
int lb_lbfluid_get_interpolated_velocity_ibm(double* p, double* v, int id);

/** Force density rescaling when external forces are defined
 */
void force_density_conversion_ibm();

#endif
#endif
