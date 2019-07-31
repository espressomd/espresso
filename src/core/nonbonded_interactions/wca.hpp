/*
  Copyright (C) 2018 The ESPResSo project

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
#ifndef WCA_HPP
#define WCA_HPP
/** \file
 *  Routines to calculate the Weeks-Chandler-Andersen potential between
 *  particle pairs.
 *
 *  Implementation in \ref wca.cpp.
 */

#include "config.hpp"

#ifdef WCA

#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <utils/math/int_pow.hpp>

int wca_set_params(int part_type_a, int part_type_b, double eps, double sig);

/** Calculate WCA force between particle p1 and p2 */
inline void add_wca_pair_force(const Particle *const p1,
                               const Particle *const p2,
                               IA_parameters *ia_params, double const d[3],
                               double dist, double force[3]) {
  if (dist < ia_params->WCA_cut) {
    auto const frac6 = Utils::int_pow<6>(ia_params->WCA_sig / dist);
    auto const fac =
        48.0 * ia_params->WCA_eps * frac6 * (frac6 - 0.5) / (dist * dist);
    for (int j = 0; j < 3; j++)
      force[j] += fac * d[j];
  }
}

/** Calculate WCA energy between particle p1 and p2. */
inline double wca_pair_energy(const Particle *p1, const Particle *p2,
                              const IA_parameters *ia_params, const double d[3],
                              double dist) {
  if (dist < ia_params->WCA_cut) {
    auto const frac6 = Utils::int_pow<6>(ia_params->WCA_sig / dist);
    return 4.0 * ia_params->WCA_eps * (Utils::sqr(frac6) - frac6 + .25);
  }
  return 0.0;
}

#endif /* ifdef WCA */
#endif
