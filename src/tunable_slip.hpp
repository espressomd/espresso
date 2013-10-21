/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2009,2010 
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
#ifndef _TUNABLE_SLIP_H
#define _TUNABLE_SLIP_H

/** \file tunable_slip.hpp
 *  Routines to generate tunable-slip boundary conditions.
 *  J.Smiatek, M.P. Allen, F. Schmid:
 *  "Tunable-slip boundaries for coarse-grained simulations of fluid flow", Europ. Phys. J. E 26, 115 (2008) 
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"

#ifdef TUNABLE_SLIP

int tunable_slip_set_params(int part_type_a, int part_type_b,
			    double temp, double gamma, double r_cut,
			    double time, double vx, double vy, double vz);

void add_tunable_slip_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params, double d[3], double dist, double force[3]);

#endif

#endif
