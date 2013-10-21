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
/** \file comforce.hpp
 *  Routines for center of mass forces
 */

#ifndef COMFORCE_H
#define COMFORCE_H

#include "utils.hpp"
#include "particle_data.hpp"
#include "statistics.hpp"

#ifdef COMFORCE

///
int comforce_set_params(int part_type_a, int part_type_b,
			int flag, int dir, double force, double fratio);

/// calculate the center of mass forces
void calc_comforce();

#endif

#endif
