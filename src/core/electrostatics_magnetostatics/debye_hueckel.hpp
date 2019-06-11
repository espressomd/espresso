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
#ifndef DEBYE_HUECKEL_H
#define DEBYE_HUECKEL_H
/** \file
 *  Routines to calculate the Debye_Hueckel  Energy or/and Debye_Hueckel force
 *  for a particle pair.
 */
#include "config.hpp"

#ifdef ELECTROSTATICS

#include "particle_data.hpp"

/** Structure to hold Debye-Hueckel Parameters. */
typedef struct {
  /** Cutoff for Debye-Hueckel interaction. */
  double r_cut;
  /** Debye kappa (inverse Debye length) . */
  double kappa;
} Debye_hueckel_params;

/** Structure containing the Debye-Hueckel parameters. */
extern Debye_hueckel_params dh_params;

/** \name Functions */
/************************************************************/
/*@{*/

int dh_set_params(double kappa, double r_cut);

/** Computes the Debye_Hueckel pair force and adds this
    force to the particle forces.
    @param q1q2      Product of the charges on p1 and p2.
    @param d         Vector pointing from p1 to p2.
    @param dist      Distance between p1 and p2.
    @param force     returns the force on particle 1.
*/
inline void add_dh_coulomb_pair_force(double const q1q2, double const d[3],
                                      double const dist, double force[3]) {
  if (dist < dh_params.r_cut) {
    double fac;
    if (dh_params.kappa > 0.0) {
      /* debye hueckel case: */
      double kappa_dist = dh_params.kappa * dist;
      fac =
          q1q2 * (exp(-kappa_dist) / (dist * dist * dist)) * (1.0 + kappa_dist);
    } else {
      /* pure Coulomb case: */
      fac = q1q2 / (dist * dist * dist);
    }
    for (int j = 0; j < 3; j++)
      force[j] += fac * d[j];
  }
}

inline double dh_coulomb_pair_energy(double const q1q2, double const dist) {
  if (dist < dh_params.r_cut) {
    if (dh_params.kappa > 0.0)
      return q1q2 * exp(-dh_params.kappa * dist) / dist;

    return q1q2 / dist;
  }
  return 0.0;
}
/*@}*/
#endif

#endif
