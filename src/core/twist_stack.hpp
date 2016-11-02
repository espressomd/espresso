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

#ifndef __TWIST_STACK_HPP
#define __TWIST_STACK_HPP

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

#ifdef TWIST_STACK
int twist_stack_set_params(int bond_type, DoubleList *params);

int calc_twist_stack_energy(Particle *si1, Particle *bi1, Particle *bi2, Particle *si2,
				      Particle *sj1, Particle *bj1, Particle *bj2, Particle *sj2,
				Bonded_ia_parameters *iaparams, double *_energy);

int calc_twist_stack_force(Particle *si1, Particle *bi1, Particle *bi2, Particle *si2,
				      Particle *sj1, Particle *bj1, Particle *bj2, Particle *sj2,
				      Bonded_ia_parameters *iaparams,
				      double f_si1[3], double f_bi1[3], double f_bi2[3], double f_si2[3],
			       double f_sj1[3], double f_bj1[3], double f_bj2[3], double f_sj2[3]);


#endif /* TWIST_STACK */

#endif /* __TWIST_STACK_HPP */
