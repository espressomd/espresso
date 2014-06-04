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
#ifndef STEPPOT_H
#define STEPPOT_H

/** \file steppot.hpp
 *  Routines to calculate the smooth step potential energy and/or force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"

#ifdef SMOOTH_STEP

///
int smooth_step_set_params(int part_type_a, int part_type_b,
			   double d, int n, double eps,
			   double k0, double sig,
			   double cut);

/** Calculate smooth step force between particle p1 and p2 */
inline void add_SmSt_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				  double d[3], double dist,double dist2, double force[3])
{
  if (!CUTOFF_CHECK(dist < ia_params->SmSt_cut))
    return;
  
  int j;
  double frac, fracP, fac=0.0,er;
  frac = ia_params->SmSt_d/dist;
  fracP = pow(frac,ia_params->SmSt_n);
  er=exp(2.*ia_params->SmSt_k0*(dist-ia_params->SmSt_sig));
  fac   =  (ia_params->SmSt_n * fracP+2.*ia_params->SmSt_eps*ia_params->SmSt_k0*dist*er/SQR(1.0+er))/dist2;

  for(j=0;j<3;j++)
    force[j] += fac * d[j];
}

/** calculate smooth step potential energy between particle p1 and p2. */
inline double SmSt_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				 double d[3], double dist,double dist2)
{
  if (!CUTOFF_CHECK(dist < ia_params->SmSt_cut))
    return 0.0;

  double frac, fracP, er;
 
  frac = ia_params->SmSt_d/dist;
  fracP = pow(frac,ia_params->SmSt_n);
  er=exp(2.*ia_params->SmSt_k0*(dist-ia_params->SmSt_sig));
  
  return fracP+ia_params->SmSt_eps/(1.0+er);
}

#endif /* ifdef SMOOTH_STEP */
#endif
