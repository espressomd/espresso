/*
  Copyright (C) 2010,2012,2013,2016 The ESPResSo project
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
#ifndef _COS2_H
#define _COS2_H

/** \file cos2.hpp
 *  Routines to calculate a flat potential with cosine tail energy and/or  force 
 *  for a particle pair.  Cosine tail is different from that in ljcos.hpp
 *  Used for attractive tail/tail interactions in lipid bilayer calculations.
 *  Same potential as ljcos2 without Lennard-Jones part.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"

#ifdef COS2
#include <cmath>

int cos2_set_params(int part_type_a, int part_type_b,
		      double eps, double offset, double w);

/** Calculate cos2 force between particle p1 and p2 */
inline void add_cos2_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  if((dist < ia_params->COS2_cut)) { 
    double fac;
    if (dist > ia_params->COS2_offset) {
      fac = -ia_params->COS2_eps*M_PI/2/ia_params->COS2_w/dist * sin(M_PI*(dist-ia_params->COS2_offset)/ia_params->COS2_w);
    }
    for(int j=0;j<3;j++)
      force[j] += fac * d[j];
  }
}

/** calculate cos2 energy between particle p1 and p2. */
inline double cos2_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{

  if((dist < ia_params->COS2_cut)) {
    if (dist < ia_params->COS2_offset)
      return -ia_params->COS2_eps;
    else
      return -ia_params->COS2_eps/2 * (cos(M_PI*(dist-ia_params->COS2_offset)/ia_params->COS2_w)+1);
  }
  return 0.0;
}

#endif /* ifdef COS2 */
#endif
