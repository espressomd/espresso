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
#ifndef _THOLE_H
#define _THOLE_H
/** \file
 *  Routines to calculate the Thole damping energy and/or force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "config.hpp"

#ifdef THOLE
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

int thole_set_params(int part_type_a, int part_type_b, double scaling_coeff,
                     double q1q2);

inline void add_thole_pair_force(const Particle *const p1,
                                 const Particle *const p2,
                                 const IA_parameters *ia_params,
                                 const double *d, double dist, double *force) {
  auto const thole_q1q2 = ia_params->THOLE_q1q2;
  auto const thole_s = ia_params->THOLE_scaling_coeff;

  if (thole_s != 0 && thole_q1q2 != 0 &&
      !(pair_bond_enum_exists_between(p1, p2, BONDED_IA_THERMALIZED_DIST))) {
    // Calc damping function (see doi.org/10.1016/0301-0104(81)85176-2)
    // S(r) = 1.0 - (1.0 + thole_s*r/2.0) * exp(-thole_s*r);
    // Calc F = - d/dr ( S(r)*q1q2/r) =
    // -(1/2)*(-2+(r^2*s^2+2*r*s+2)*exp(-s*r))*q1q2/r^2 Everything before
    // q1q2/r^2 can be used as a factor for the Coulomb::central_force method
    double sr = thole_s * dist;
    double dS_r = 0.5 * (2.0 - (exp(-sr) * (sr * (sr + 2.0) + 2.0)));
    auto const f = Coulomb::central_force(thole_q1q2 * (-1. + dS_r), d, dist);

    force[0] += f[0];
    force[1] += f[1];
    force[2] += f[2];
  }
}

inline double thole_pair_energy(const Particle *p1, const Particle *p2,
                                const IA_parameters *ia_params,
                                const double d[3], double dist) {
  double thole_s = ia_params->THOLE_scaling_coeff;
  double thole_q1q2 = ia_params->THOLE_q1q2;

  if (thole_s != 0 && thole_q1q2 != 0 && dist < Coulomb::cutoff(box_l) &&
      !(pair_bond_enum_exists_between(p1, p2, BONDED_IA_THERMALIZED_DIST))) {

    // Add damped p3m shortrange energy
    double sd = thole_s * dist;
    double S_r = 1.0 - (1.0 + sd / 2.0) * exp(-sd);
    return Coulomb::pair_energy(p1, p2, thole_q1q2 * (-1. + S_r), d,
                                dist * dist, dist);
  }
  return 0.0;
}
#endif
#endif
