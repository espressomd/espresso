/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

/** \file constraint.hpp
 *  Routines for handling of constraints.
 *  Only active if the feature CONSTRAINTS is activated.
 *  see also \ref interaction_data.hpp
 */

#include "particle_data.hpp"
#include "interaction_data.hpp"

#ifdef CONSTRAINTS

extern int reflection_happened;

/** Exported functions
 */

Constraint *generate_constraint();


void reflect_particle(Particle *p1, double *distance_vec, 
		      int reflecting);

void add_constraints_forces(Particle *p1);

double add_constraints_energy(Particle *p1);

void init_constraint_forces();
#endif

#endif
