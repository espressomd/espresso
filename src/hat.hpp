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
#ifndef hat_H
#define hat_H

/** \file soft_sphere.hpp
 *  Routines to calculate the soft-sphere energy and/or  force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"

#ifdef HAT

///
int hat_set_params(int part_type_a, int part_type_b,
			   double Fmax, double r);

/** Resultant Force due to an hat potential between two
    particles at interatomic separation dist */
inline double hat_force_r(double Fmax, double r, double dist )
{
  return dist < r ? Fmax * (1 - dist/r) : 0.0;
}

/** Potential Energy due to an hat potential between two
    particles at interatomic separation dist */
inline double hat_energy_r(double Fmax, double r, double dist )
{
  return dist < r ? Fmax * (dist-r) * ( (dist+r) / (2.0*r) - 1.0 ) : 0.0;
}

/** Calculate hat potential force between particle p1 and p2 */
inline void add_hat_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double fac=0.0;
  if(CUTOFF_CHECK(dist < ia_params->HAT_r)) {
    if(dist < ia_params->HAT_r) {
      fac = hat_force_r(ia_params->HAT_Fmax, ia_params->HAT_r, dist)/dist;
      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: Hat-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    
    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: hat    f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: hat    f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

/** calculate hat energy between particle p1 and p2. */
inline double hat_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  if(CUTOFF_CHECK(dist < ia_params->HAT_r)) {   
    return hat_energy_r(ia_params->HAT_Fmax, ia_params->HAT_r, dist);    
  }
  return 0.0;
}

#endif /* ifdef HAT */
#endif
