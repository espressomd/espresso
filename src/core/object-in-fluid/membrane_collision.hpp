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
#ifndef MEMBRANE_COLLISION_H
#define MEMBRANE_COLLISION_H

/** \file membrane_collision.hpp
 *  Routines to calculate the membrane collision force
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "../utils.hpp"
#include "../interaction_data.hpp"
#include "../particle_data.hpp"
#include "../mol_cut.hpp"
#include "../integrate.hpp"
#include "../random.hpp"
#include "../grid.hpp"

#ifdef MEMBRANE_COLLISION

int membrane_collision_set_params(int part_type_a, int part_type_b,
			   double a, double n, double cut, double offset);

//** Resultant force due to a sigmoid potential between two
// particles at interatomic separation r */
inline double sigmoid_force_r(double a, double n, double r )
{
    return (a/(1+exp(n*r)));
}

/** Calculate membrane-collision force between particle p1 and p2 */
inline void add_membrane_collision_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
	/************************
	*  
	* Description of implementation:
    * We have two particles, each belongs to the membrane of a different immersed object
	* For both particles we have the position of the particle - p, and in part->p.out_direction are the coordinates of the outward normal vector (with respect to the immersed object).
    *
	* Algorithm:
	* 1) compute angle between outward normal vectors out1 and out2 to check if they are almost collinear (in this case nothing happens, this is at the edge of contact zone and forces should not be applied because they would quite likely have wrong direction)
    * 2) in other cases, repulsive forces are applied in the direction out1-out2
    *
    *********************/
	
    int j;
    double r_off, fac=0.0, product, angle, ndir;
    double out1[3],out2[3],dir[3];
    
    if(CUTOFF_CHECK(dist < ia_params->membrane_cut+ia_params->membrane_offset)) {
        
        r_off = dist - ia_params->membrane_offset;
        // offset needs to be checked for the unphysical case when r_off should be negative
        
        if(r_off > 0.0) {
            
            memmove(out1, p1->p.out_direction, 3*sizeof(double));
            memmove(out2, p2->p.out_direction, 3*sizeof(double));
            
            // this is the direction in which the repulsive forces will be applied and its norm
            vector_subt(dir,out1,out2);
            ndir=normr(dir);
            
            // for very small angles the force should not be applied - these happen at the crossing of the boundary and would result in oscillation
            product = scalar(out1,out2);
            angle = acos(product);
            
            if (fabs(angle)>SMALL_OIF_MEMBRANE_CUTOFF){
                fac = sigmoid_force_r(ia_params->membrane_a, ia_params->membrane_n, r_off)/dist;
                for(j=0;j<3;j++)
                    force[j] -= fac * dir[j]/ndir;
            }
        }
    }
}



#endif /* ifdef MEMBRANE_COLLISION */
#endif
