// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#ifndef VIRTUAL_SITES_H
#define VIRTUAL_SITES_H

#include "particle_data.h"
#include <tcl.h>

/** \file virtual_sites.h
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
MDINLINE int ifParticleIsVirtual(Particle *p){
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
 #include "virtual_sites_com.h"
#endif

// Virtual particles relative to position and orientation of a real particle
#ifdef VIRTUAL_SITES_RELATIVE
 #include "virtual_sites_relative.h"
#endif

#endif

#endif
