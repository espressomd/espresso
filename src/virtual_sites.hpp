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
#ifndef _VIRTUAL_SITES_H
#define _VIRTUAL_SITES_H

#include "particle_data.hpp"

/** \file virtual_sites.hpp
 *  This file contains routine to handle virtual sites
 *  Virtual sites are like particles, but they will be not integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move no-virtual particles
 *  - update virtual sites
 */

#ifdef VIRTUAL_SITES
// Recalculate position and velocity for all virtual particles
void update_mol_vel_pos();
// Recalc velocities for virtual particles
void update_mol_vel();
// Recalc positions of virtual particles
void update_mol_pos();


// Update the position of all virutal particles 
// in the partCfg-array rather than in the local cells.
int update_mol_pos_cfg();


// The following three functions have to be provided by all implementations
// of virtual sites
// Update the vel/pos of the given virtual particle as defined by the real 
// particles in the same molecule
// void update_mol_pos_particle(Particle *);
// void update_mol_vel_particle(Particle *);

// Distribute forces that have accumulated on virtual particles to the 
// associated real particles
//void distribute_mol_force();


// Checks, if a particle is virtual
inline int ifParticleIsVirtual(Particle *p){
   if (p->p.isVirtual == 0) {
      return 0;
   }
   else{
      return 1;
   }
}

// According to what rules the virtual particles are placed and the forces and
// torques accumulating on the virtual particles distributed back to real 
// particles, is decided by a specific implementation.

// Virtual particles in center of mass of molecule
#ifdef VIRTUAL_SITES_COM
 #include "virtual_sites_com.hpp"
#endif

// Virtual particles relative to position and orientation of a real particle
#ifdef VIRTUAL_SITES_RELATIVE
 #include "virtual_sites_relative.hpp"
#endif

#endif

#endif
