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

#ifndef _VIRTUAL_SITES_IBM_H
#define _VIRTUAL_SITES_IBM_H

#include "config.hpp"
#include "particle_data.hpp"

extern int integration_rule_ibm;

#ifdef VIRTUAL_SITES_IMMERSED_BOUNDARY

//Update Position ~ Euler/Runge-Kutta/Adams-Bashforth
void update_mol_pos_particle(Particle *);
//Update Velocity ~ Get interpolated velocity of LB , update old velocity for ab integration
void update_mol_vel_particle(Particle *);
//Since no 'real' particles are involved, this function will stay empty
void distribute_mol_force();

#endif 

#endif
