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
#ifndef _MORSE_H
#define _MORSE_H

/** \file
 *  Routines to calculate the Morse potential between particle pairs.
 *
 *  Implementation in \ref morse.cpp.
 */

#include "config.hpp"

#ifdef MORSE

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

int morse_set_params(int part_type_a, int part_type_b, double eps, double alpha,
                     double rmin, double cut);

/** Calculate Morse force between particle p1 and p2 */
inline void add_morse_pair_force(Particle const *const p1,
                                 Particle const *const p2,
                                 IA_parameters const *const ia_params,
                                 Utils::Vector3d const &d, double dist,
                                 Utils::Vector3d &force) {
  if (dist < ia_params->morse.cut) {
    auto const add =
        exp(-ia_params->morse.alpha * (dist - ia_params->morse.rmin));
    double fac = -ia_params->morse.eps * 2.0 * ia_params->morse.alpha *
                 (add - Utils::sqr(add)) / dist;
    force += fac * d;

#ifdef MORSE_WARN_WHEN_CLOSE
    if (fac * dist > 1000)
      fprintf(stderr, "%d: Morse-Warning: Pair (%d-%d) force=%f dist=%f\n",
              this_node, p1->p.identity, p2->p.identity, fac * dist, dist);
#endif
    ONEPART_TRACE(if (p1->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: MORSE   f = (%.3e,%.3e,%.3e) "
                              "with part id=%d at dist %f fac %.3e\n",
                              this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                              p2->p.identity, dist, fac));
    ONEPART_TRACE(if (p2->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: MORSE   f = (%.3e,%.3e,%.3e) "
                              "with part id=%d at dist %f fac %.3e\n",
                              this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                              p1->p.identity, dist, fac));

    MORSE_TRACE(fprintf(
        stderr, "%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
        this_node, p1->p.identity, p2->p.identity, dist, force[0], force[1],
        force[2]));
  }
}

/** Calculate Morse energy between particle p1 and p2. */
inline double morse_pair_energy(Particle const *const p1,
                                Particle const *const p2,
                                IA_parameters const *const ia_params,
                                Utils::Vector3d const &d, double dist) {
  if (dist < ia_params->morse.cut) {
    auto const add =
        exp(-ia_params->morse.alpha * (dist - ia_params->morse.rmin));
    auto const fac = ia_params->morse.eps * (Utils::sqr(add) - 2 * add) -
                     ia_params->morse.rest;
    return fac;
  }
  return 0.0;
}

#endif /* ifdef MORSE */
#endif /* ifdef _MORSE_H */
