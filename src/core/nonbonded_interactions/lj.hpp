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
#ifndef _LJ_H
#define _LJ_H

#include "config.hpp"

#ifdef LENNARD_JONES

/** \file
 *  Routines to calculate the Lennard-Jones energy and/or  force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

int lennard_jones_set_params(int part_type_a, int part_type_b, double eps,
                             double sig, double cut, double shift,
                             double offset, double min);

/** Calculate Lennard-Jones force between particle p1 and p2 */
inline void add_lj_pair_force(const Particle *const p1,
                              const Particle *const p2,
                              IA_parameters *ia_params, double const d[3],
                              double dist, double force[3]) {
  if ((dist < ia_params->LJ_cut + ia_params->LJ_offset) &&
      (dist > ia_params->LJ_min + ia_params->LJ_offset)) {
    double r_off = dist - ia_params->LJ_offset;
    double frac2 = Utils::sqr(ia_params->LJ_sig / r_off);
    double frac6 = frac2 * frac2 * frac2;
    double fac =
        48.0 * ia_params->LJ_eps * frac6 * (frac6 - 0.5) / (r_off * dist);
    for (int j = 0; j < 3; j++)
      force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
    if (fac * dist > 1000)
      fprintf(stderr, "%d: LJ-Warning: Pair (%d-%d) force=%f dist=%f\n",
              this_node, p1->p.identity, p2->p.identity, fac * dist, dist);
#endif
    ONEPART_TRACE(if (p1->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with "
                              "part id=%d at dist %f fac %.3e\n",
                              this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                              p2->p.identity, dist, fac));
    ONEPART_TRACE(if (p2->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with "
                              "part id=%d at dist %f fac %.3e\n",
                              this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                              p1->p.identity, dist, fac));

    LJ_TRACE(fprintf(
        stderr, "%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
        this_node, p1->p.identity, p2->p.identity, dist, fac * d[0], fac * d[1],
        fac * d[2]));
  }
}

/** calculate Lennard-Jones energy between particle p1 and p2. */
inline double lj_pair_energy(const Particle *p1, const Particle *p2,
                             const IA_parameters *ia_params, const double d[3],
                             double dist) {
  if ((dist < ia_params->LJ_cut + ia_params->LJ_offset) &&
      (dist > ia_params->LJ_min + ia_params->LJ_offset)) {
    double r_off = dist - ia_params->LJ_offset;
    double frac2 = Utils::sqr(ia_params->LJ_sig / r_off);
    double frac6 = frac2 * frac2 * frac2;
    return 4.0 * ia_params->LJ_eps *
           (Utils::sqr(frac6) - frac6 + ia_params->LJ_shift);
  }
  return 0.0;
}

#endif /* ifdef LENNARD_JONES */
#endif
