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


#include <list>
#include <string>

#include "particle_data.hpp"

namespace Scafacos {
#if defined(SCAFACOS) 
/** Near-field pair force */
void add_pair_force(Particle *p1, Particle *p2, double *d, double dist, double *force);
/** Near-field pair energy */
double pair_energy(Particle* p1, Particle* p2, double dist);
/** Long range part */
void add_long_range_force();
/** Calculate long range energy contribution */
double long_range_energy(); 
/** Get parameters */
std::string get_method_and_parameters();
/** Set parameters */
void set_parameters(const std::string &method, const std::string &params, bool dipolar);
double get_r_cut();

/** Is scafacos used for dipolar interactions */
bool dipolar();

/** Choose whether scafacos is used for dipolar interactions */
void set_dipolar(bool d);

/** Reinit scafacos number of particles, box shape and periodicity */
void update_system_params();

#endif /* SCAFACOS */

std::list<std::string> available_methods();

void free_handle();

} // namespace scafacos

/** Parameter callback */
void mpi_scafacos_set_parameters_slave(int n_method, int n_params);
void mpi_scafacos_set_r_cut_and_tune_slave(int a, int b);
void mpi_scafacos_free_slave(int a, int b) ;

#endif
