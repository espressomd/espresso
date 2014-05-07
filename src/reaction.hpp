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
#ifndef REACTION_H
#define REACTION_H
/** \file reaction.hpp
 *
 */
 
#include "utils.hpp"
#include "particle_data.hpp"

typedef struct {
	int reactant_type;
	int product_type;
	int catalyzer_type;
	double range;
	double ct_rate;
	double eq_rate;
  int sing_mult;
}  reaction_struct;

extern reaction_struct reaction;

#ifdef CATALYTIC_REACTIONS
/** sanity checks for the reaction code */
void reactions_sanity_checks();
/** broadcasts reaction parameters and sets up an entry in the ia_params, so
    that the verlet radius is equal or bigger than the reaction range.
**/
void local_setup_reaction();
void integrate_reaction();
#endif

#endif /* ifdef REACTION_H */
