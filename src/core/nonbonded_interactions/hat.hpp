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
#ifndef hat_H
#define hat_H

/** \file
 *  Routines to calculate the soft-sphere energy and/or  force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "config.hpp"

#ifdef HAT

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

///
int hat_set_params(int part_type_a, int part_type_b, double Fmax, double r);

/** Resultant Force due to an hat potential between two
    particles at interatomic separation dist */
inline double hat_force_r(double Fmax, double r, double dist) {
  return dist < r ? Fmax * (1 - dist / r) : 0.0;
}

/** Potential Energy due to an hat potential between two
    particles at interatomic separation dist */
inline double hat_energy_r(double Fmax, double r, double dist) {
  return dist < r ? Fmax * (dist - r) * ((dist + r) / (2.0 * r) - 1.0) : 0.0;
}

/** Calculate hat potential force between particle p1 and p2 */
inline void add_hat_pair_force(const Particle *const p1,
                               const Particle *const p2,
                               IA_parameters *ia_params, double const d[3],
                               double dist, double force[3]) {
  if (dist > 0. && dist < ia_params->HAT_r) {
    auto const fac =
        hat_force_r(ia_params->HAT_Fmax, ia_params->HAT_r, dist) / dist;
    for (int j = 0; j < 3; j++)
      force[j] += fac * d[j];
  }
}

/** calculate hat energy between particle p1 and p2. */
inline double hat_pair_energy(const Particle *p1, const Particle *p2,
                              const IA_parameters *ia_params, const double d[3],
                              double dist) {
  if ((dist < ia_params->HAT_r)) {
    return hat_energy_r(ia_params->HAT_Fmax, ia_params->HAT_r, dist);
  }
  return 0.0;
}

#endif /* ifdef HAT */
#endif
