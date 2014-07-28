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

#include "immersed-boundary/lbtracers.hpp"
#include "lb.hpp"
#include "integrate.hpp"
#include "initialize.hpp"
#include "communication.hpp"
#include "particle_data.hpp"

#ifdef LBTRACERS

//Update Position ~ Euler/Adams Bashforth
void update_mol_pos_particle(Particle *p) {
  //Do Euler/Adams Bashforth for particle p; assume velocity and previous velocity has already been calculated 
  // & is stored in particle data
  int j;
  double skin2 = SQR(0.5 * skin);
		
  for(j=0;j<3;j++) 
    {
      // // Euler
      // p->r.p[j] = p->r.p[j] + p->m.v[j]*time_step;

      // Two-step Adams–Bashforth
      p->r.p[j] = p->r.p[j] + (((1.5*p->m.v[j]) - (0.5*p->l.v_old[j]))*time_step);
    }
	    
  // Check if a particle might have crossed a box border (Verlet criterium); 
  //if possible resort_particles = 1
  const double dist2 = distance2( p->r.p, p->l.p_old);	
  if ( dist2 > skin2 ) { resort_particles = 1; }		
}

//Update Velocity ~ Get interpolated velocity of LB set old velocity of particle
void update_mol_vel_particle(Particle *p) 
{
  int j;
  double p_temp[3];
	
  for(j=0;j<3;j++) {
    p_temp[j] = p->r.p[j];
  }
	
  if ( !(lattice_switch & LATTICE_LB_GPU) )
    {
      // Need to interpolate velocity here only for CPU
      // For GPU it is already stored
      double v_int[3] = {0,0,0};
      lb_lbfluid_get_interpolated_velocity_lbtrace(p_temp,v_int, p->p.identity);
            
      for ( j = 0; j < 3; j++){ 
       
	p->l.v_old[j] = p->m.v[j];
	p->m.v[j] = v_int[j];

      }
    }
}

//Distribute forces
void distribute_mol_force() {
  //All forces these sites are subject to are influencing the LB fluid, not other
  //particles, therefore no forces need to be distributed here. => Do nothing
}


#endif
