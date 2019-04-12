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
#ifndef REACTION_FIELD_H
#define REACTION_FIELD_H
/** \file
 *  Routines to calculate the Reaction Field Energy or/and force
 *  for a particle pair.
 *  M. Neumann, J. Chem. Phys 82, 5663 (1985)
 *  \ref forces.cpp
 *
 */

#include "config.hpp"

#ifdef ELECTROSTATICS
#include "particle_data.hpp"

/** Structure to hold Reaction Field Parameters. */
typedef struct {
  /** ionic strength . */
  double kappa;
  /** epsilon1 (continuum dielectric constant inside) . */
  double epsilon1;
  /** epsilon2 (continuum dielectric constant outside) . */
  double epsilon2;
  /** Cutoff for Reaction Field interaction. */
  double r_cut;
  /** B important prefactor . */
  double B;
} Reaction_field_params;

/** Structure containing the Reaction Field parameters. */
extern Reaction_field_params rf_params;

/** \name Functions */
/************************************************************/
/*@{*/

///
int rf_set_params(double kappa, double epsilon1, double epsilon2, double r_cut);

inline void add_rf_coulomb_pair_force_no_cutoff(double const q1q2,
                                                double const d[3],
                                                double const dist,
                                                double force[3]) {
  int j;
  double fac;
  fac = 1.0 / (dist * dist * dist) +
        rf_params.B / (rf_params.r_cut * rf_params.r_cut * rf_params.r_cut);
  fac *= q1q2;

  for (j = 0; j < 3; j++)
    force[j] += fac * d[j];
}

/** Computes the Reaction Field pair force and adds this
    force to the particle forces.
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param d         Vector pointing from p1 to p2.
    @param dist      Distance between p1 and p2.
    @param force     returns the force on particle 1.
*/
inline void add_rf_coulomb_pair_force(double const q1q2, double const d[3],
                                      double const dist, double force[3]) {
  if (dist < rf_params.r_cut) {
    add_rf_coulomb_pair_force_no_cutoff(q1q2, d, dist, force);
  }
}

inline double rf_coulomb_pair_energy_no_cutoff(double const q1q2,
                                               double const dist) {
  double fac;
  fac = 1.0 / dist -
        (rf_params.B * dist * dist) /
            (2 * rf_params.r_cut * rf_params.r_cut * rf_params.r_cut);
  // cut off part
  fac -= (1 - rf_params.B / 2) / rf_params.r_cut;
  fac *= q1q2;
  return fac;
}

inline double rf_coulomb_pair_energy(double const q1q2, double const dist) {
  if (dist < rf_params.r_cut) {
    return rf_coulomb_pair_energy_no_cutoff(q1q2, dist);
  }
  return 0.0;
}

/*@}*/
#endif

#endif
