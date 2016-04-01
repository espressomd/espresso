/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

/** \file This file contains the c-type wrapper interface to the (oop-) scafacos interface. */

#ifndef __SCAFACOS_HPP
#define __SCAFACOS_HPP

#include "config.hpp"

#if defined(SCAFACOS) and defined(ELECTROSTATICS)

#include <list>
#include <string>

#include "particle_data.hpp"
#endif /* SCAFACOS */

namespace Electrostatics {
namespace Scafacos {
#if defined(SCAFACOS) and defined(ELECTROSTATICS)
/** Pair force */
void add_pair_force(Particle *p1, Particle *p2, double *d, double dist, double *force);
/** Long range part */
void add_long_range_force();
/** Get parameters */
std::string get_parameters();
/** Set parameters */
void set_parameters(const std::string &method, const std::string &params);
double get_r_cut();
std::list<std::string> available_methods();
#endif /* SCAFACOS */

}
}

/** Parameter callback */
void mpi_scafacos_set_parameters_slave(int n_method, int n_params);

#endif
