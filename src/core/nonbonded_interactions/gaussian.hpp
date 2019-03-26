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
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

/** \file
 *  Routines to calculate the Gaussian energy and/or force
 *  for a particle pair.
 */

#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

#ifdef GAUSSIAN

///
int gaussian_set_params(int part_type_a, int part_type_b, double eps,
                        double sig, double cut);

/** Calculate Gaussian force between particle p1 and p2 */
inline void add_gaussian_pair_force(const Particle *const p1,
                                    const Particle *const p2,
                                    IA_parameters *ia_params, double const d[3],
                                    double dist, double dist2,
                                    double force[3]) {
  double fac;
  int j;
  if ((dist < ia_params->Gaussian_cut)) {
    fac = ia_params->Gaussian_eps / pow(ia_params->Gaussian_sig, 2) *
          exp(-0.5 * Utils::sqr(dist / ia_params->Gaussian_sig));

    for (j = 0; j < 3; j++)
      force[j] += fac * d[j];
  }
}

/** calculate Lennard-Jones energy between particle p1 and p2. */
inline double gaussian_pair_energy(const Particle *p1, const Particle *p2,
                                   const IA_parameters *ia_params,
                                   const double d[3], double dist,
                                   double dist2) {
  if ((dist < ia_params->Gaussian_cut)) {
    return ia_params->Gaussian_eps *
           exp(-0.5 * Utils::sqr(dist / ia_params->Gaussian_sig));
  }
  return 0.0;
}

#endif /* ifdef GAUSSIAN */
#endif
