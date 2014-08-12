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

#ifndef CG_DNA_HPP
#define CG_DNA_HPP

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

#ifdef CG_DNA

// extrema for cos(theta), used for the force calculations that involve angles
#define COS_MAX (0.99999999)
#define COS_MIN (-0.99999999)

int cg_dna_basepair_set_params(int bond_type, DoubleList *params);
int cg_dna_stacking_set_params(int bond_type, DoubleList *params);

int calc_cg_dna_stacking_energy(Particle *si1, Particle *bi1, Particle *bi2, Particle *si2,
				      Particle *sj1, Particle *bj1, Particle *bj2, Particle *sj2,
				Bonded_ia_parameters *iaparams, double *_energy);

int calc_cg_dna_stacking_force(Particle *si1, Particle *bi1, Particle *bi2, Particle *si2,
				      Particle *sj1, Particle *bj1, Particle *bj2, Particle *sj2,
				      Bonded_ia_parameters *iaparams,
				      double f_si1[3], double f_bi1[3], double f_bi2[3], double f_si2[3],
			       double f_sj1[3], double f_bj1[3], double f_bj2[3], double f_sj2[3]);

int calc_cg_dna_basepair_force(Particle *s1, Particle *b1, Particle *b2, Particle *s2, Bonded_ia_parameters *iaparams, double f_s1[3], double f_b1[3], double f_b2[3], double f_s2[3]);

int calc_cg_dna_basepair_energy(Particle *s1, Particle *b1, Particle *b2, Particle *s2, Bonded_ia_parameters *iaparams, double *_energy);

#endif /* CG_DNA */

#endif /* CG_DNA_HPP */
