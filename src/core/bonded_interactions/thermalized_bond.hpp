
/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#ifndef THERMALIZED_DIST_H
#define THERMALIZED_DIST_H
/** \file
 *  Routines to thermalize the com and distance of a particle pair.
 *  \ref forces.cpp
 */

/** number of thermalized bonds */
extern int n_thermalized_bonds;

#include "bonded_interaction_data.hpp"
#include "debug.hpp"
#include "integrate.hpp"
#include "random.hpp"
#include "utils.hpp"

// Set the parameters for the thermalized bond
int thermalized_bond_set_params(int bond_type, double temp, double gamma, double r_cut);

void thermalized_bond_heat_up();
void thermalized_bond_cool_down();
void thermalized_bond_update_params(double pref_scale);
void thermalized_bond_init();

/** Separately thermalizes the com and distance of a particle pair
    and adds this force to the particle forces.
    @param p1        First particle
    @param p2        Second/middle particle
    @param iaparams  Parameters of interaction
    @param dx        Change in position
    @param force1    Force on particle p1
    @param force2    Force on particle p2
    @return true if bond is broken
*/

inline int calc_thermalized_bond_forces(const Particle *p1, const Particle *p2,
                                        const Bonded_ia_parameters *iaparams,
                                        double dx[3], double force1[3],
                                        double force2[3]) {
  // Bond broke?
  if (iaparams->p.thermalized_bond.r_cut > 0.0 &&
      Vector3d(dx, dx + 3).norm() > iaparams->p.thermalized_bond.r_cut) {
    return 1;
  }

  double force_lv_dist, dist_vel;
  double mass_tot = p1->p.mass + p2->p.mass;
  double sqrt_mass_red = sqrt(p1->p.mass * p2->p.mass / mass_tot);

  for (int i = 0; i < 3; i++) {
    // Langevin thermostat for distance p1->p2
    dist_vel = p2->m.v[i] - p1->m.v[i];
    if (iaparams->p.thermalized_bond.pref2 > 0.0) {
      force_lv_dist = -iaparams->p.thermalized_bond.pref1 * dist_vel +
                      sqrt_mass_red * iaparams->p.thermalized_bond.pref2 *
                          (d_random() - 0.5);
    } else {
      force_lv_dist = -iaparams->p.thermalized_bond.pref1 * dist_vel;
    }
    // Add forces
    force1[i] = - force_lv_dist;
    force2[i] = + force_lv_dist;
  }
  return 0;
}

#endif
