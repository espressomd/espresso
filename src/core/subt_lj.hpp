/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef SUBT_LJ_H
#define SUBT_LJ_H
/** \file subt_lj.hpp
 *  Routines to subtract the LENNARD-JONES Energy and/or the LENNARD-JONES force
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "config.hpp"

#ifdef LENNARD_JONES

#include "debug.hpp"
#include "interaction_data.hpp"
#include "lj.hpp"
#include "utils.hpp"

/// set the parameters for the subtract LJ potential
int subt_lj_set_params(int bond_type);

/** Computes the negative of the LENNARD-JONES pair forces
    and adds this force to the particle forces.
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  Parameters of interaction
    @param dx        change in position
    @param force     force on particles
    @return true if bond is broken
*/
inline int calc_subt_lj_pair_force(Particle *p1, Particle *p2,
                                   Bonded_ia_parameters *iaparams, double dx[3],
                                   double force[3]) {
  auto ia_params = get_ia_param(p1->p.type, p2->p.type);

  for (int i = 0; i < 3; i++) {
    dx[i] *= -1;
  }

  add_lj_pair_force(p1, p2, ia_params, dx, Utils::veclen(dx), force);

  return ES_OK;
}

inline int subt_lj_pair_energy(Particle *p1, Particle *p2,
                               Bonded_ia_parameters *iaparams, double dx[3],
                               double *_energy) {
  auto ia_params = get_ia_param(p1->p.type, p2->p.type);

  *_energy = -lj_pair_energy(p1, p2, ia_params, dx, Utils::veclen(dx));
  return ES_OK;
}

#endif

#endif
