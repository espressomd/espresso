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
#ifndef BMHTF_NACL_H
#define BMHTF_NACL_H
/** \file bmhtf-nacl.hpp
 *  Routines to calculate the Born-Meyer-Huggins-Tosi-Fumi energy and/or force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"

#ifdef BMHTF_NACL

///
int BMHTF_set_params(int part_type_a, int part_type_b,
		     double A, double B, double C,
		     double D, double sig, double cut);

/** Calculate smooth step force between particle p1 and p2 */
inline void add_BMHTF_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				   double d[3], double dist, double dist2, double force[3])
{
  int j;
  double pw8, fac = 0.0;
  if(CUTOFF_CHECK(dist < ia_params->BMHTF_cut)) {
    pw8 = dist2*dist2*dist2*dist2;
    fac = ia_params->BMHTF_A*ia_params->BMHTF_B*
      exp(ia_params->BMHTF_B*(ia_params->BMHTF_sig - dist))/dist -
      6*ia_params->BMHTF_C/pw8 - 8*ia_params->BMHTF_D/pw8/dist2;

    for(j=0;j<3;j++) force[j] += fac * d[j];
  }
}

/** calculate smooth step potential energy between particle p1 and p2. */
inline double BMHTF_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				  double d[3], double dist, double dist2)
{
  double pw6;
 
  if(CUTOFF_CHECK(dist < ia_params->BMHTF_cut)) {
    pw6 = dist2*dist2*dist2;
    return ia_params->BMHTF_A*
      exp(ia_params->BMHTF_B*(ia_params->BMHTF_sig - dist)) -
      ia_params->BMHTF_C/pw6 - ia_params->BMHTF_D/pw6/dist2 + ia_params->BMHTF_computed_shift;
  }
  return 0.0;
}

#endif
#endif
