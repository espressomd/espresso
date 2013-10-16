/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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
#ifndef _OBJECT_IN_FLUID_STRETCHING_FORCE_H
#define _OBJECT_IN_FLUID_STRETCHING_FORCE_H

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"

/** \file stretching_force.hpp
 *  Routines to calculate the STRETCHING_FORCE Energy or/and STRETCHING_FORCE force 
 *  for a particle pair. (Dupin2007)
 *  \ref forces.cpp
*/

/************************************************************/

/// set the parameters for the stretching_force potential
int stretching_force_set_params(int bond_type, double r0, double ks);


inline double KS(double lambda){ // Defined by (19) from Dupin2007
	double res;
	res = (pow(lambda,0.5) + pow(lambda,-2.5))/(lambda + pow(lambda,-3.));
	return res;
}

/** Computes the STRETCHING_FORCE pair force and adds this
 force to the particle forces (see \ref tclcommand_inter).
 @param p1        Pointer to first particle.
 @param p2        Pointer to second particle.
 @param iaparams  spring stiffness ks, initial distance between particles (see \ref tclcommand_inter).
 @param dx        particle distance vector
 @param force     returns force of particle 1
 @return true if the bond is broken
 */
inline int calc_stretching_force_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac, dr, len2, len;
  double lambda;

  len2 = sqrlen(dx);

  len = sqrt(len2);
  dr = len - iaparams->p.stretching_force.r0;

  lambda = 1.0*len/iaparams->p.stretching_force.r0;
  fac = -iaparams->p.stretching_force.ks * KS(lambda) * dr / iaparams->p.stretching_force.r0;
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i]/len;

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(len2),fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(len2),fac));
  
  return 0;
}

#endif
