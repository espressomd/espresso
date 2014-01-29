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
#ifndef umbrella_H
#define umbrella_H

/** \file umbrella.hpp
 *  Routines to calculate the umbrella energy and/or  force 
 *  for a particle pair: harmonic interaction in only one 
 *  cartesian direction. Useful for umbrella sampling simulations.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"

#ifdef UMBRELLA

///
int umbrella_bonded_set_params(int bond_type, double k,
			   int dir, double r);


/** Resultant Force due to an umbrella potential */
inline double umbrella_force_r(double k, int dir, double r, double distn )
{
  return -k * (distn - r);
}

/** Potential Energy due to an umbrella potential */
inline double umbrella_energy_r(double k, int dir, double r, double distn )
{
  return 0.5 * k * SQR(distn - r);
}

/** Calculate umbrella potential force between particle p1 and p2 */
inline void add_umbrella_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  double distn;
  double fac=0.0;
  distn = d[ia_params->UMBRELLA_dir];
  fac = umbrella_force_r(ia_params->UMBRELLA_k, ia_params->UMBRELLA_dir,
			 ia_params->UMBRELLA_r, distn);
  force[ia_params->UMBRELLA_dir] += fac;
      
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: umbrella f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: umbrella f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  
}

/** calculate umbrella energy between particle p1 and p2. */
inline double umbrella_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{ 
  double distn;
  distn = d[ia_params->UMBRELLA_dir];  
  return umbrella_energy_r(ia_params->UMBRELLA_k, ia_params->UMBRELLA_dir,
			   ia_params->UMBRELLA_r, distn);    
}

#endif /* ifdef UMBRELLA */
#endif
