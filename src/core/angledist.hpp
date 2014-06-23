/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#ifndef _ANGLEDIST_H
#define _ANGLEDIST_H
/** \file angledist.hpp
 *  Routines to calculate the angle and distance dependent (from a constraint) energy or/and and force
 *  for a particle triple.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"

#ifdef BOND_ANGLEDIST

/** set parameters for the angledist potential. The type of the
    angledist potential is chosen via config.hpp and cannot be changed
    at runtime.
**/
int angledist_set_params(int bond_type, double bend, double phimin,
			 double distmin, double phimax, double distmax);

///
int calc_angledist_force(Particle *p_mid, Particle *p_left,
			 Particle *p_right,
			 Bonded_ia_parameters *iaparams,
			 double force1[3], double force2[3]);

/** Computes the three body angle interaction energy (see \ref tclcommand_inter, \ref tclcommand_analyze). 
    @param p_mid        Pointer to second/middle particle.
    @param p_left       Pointer to first particle.
    @param p_right      Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param _energy   return energy pointer.
    @return 0.
*/
int angledist_energy(Particle *p_mid, Particle *p_left, Particle *p_right, 
		     Bonded_ia_parameters *iaparams, double *_energy);

#endif /* BOND_ANGLEDIST */
#endif /* ANGLEDIST_H */
