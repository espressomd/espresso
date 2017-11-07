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
#ifndef _LJCOS2_H
#define _LJCOS2_H

/** \file ljcos2.hpp
 *  Routines to calculate the lennard-jones with cosine tail energy and/or force
 *  for a particle pair.  Cosine tail is different from that in ljcos.hpp
 *  Used for attractive tail/tail interactions in lipid bilayer calculations
 *  \ref forces.cpp
*/

#include "config.hpp"

#ifdef LJCOS2

#include "debug.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

#include <cmath>

int ljcos2_set_params(int part_type_a, int part_type_b, double eps, double sig,
                      double offset, double w);

/** Calculate lj-cos2 force between particle p1 and p2 */
inline void add_ljcos2_pair_force(const Particle *const p1,
                                  const Particle *const p2,
                                  IA_parameters *ia_params, double d[3],
                                  double dist, double force[3]) {
  int j;
  double r_off, frac2, frac6, fac = 0.0;
  if ((dist < ia_params->LJCOS2_cut + ia_params->LJCOS2_offset)) {
    r_off = dist - ia_params->LJCOS2_offset;
    if (r_off < ia_params->LJCOS2_rchange) {
      frac2 = SQR(ia_params->LJCOS2_sig / r_off);
      frac6 = frac2 * frac2 * frac2;
      fac =
          48.0 * ia_params->LJCOS2_eps * frac6 * (frac6 - 0.5) / (r_off * dist);
    } else if (r_off < ia_params->LJCOS2_rchange + ia_params->LJCOS2_w) {
      fac =
          -ia_params->LJCOS2_eps * M_PI / 2 / ia_params->LJCOS2_w / dist *
          sin(M_PI * (r_off - ia_params->LJCOS2_rchange) / ia_params->LJCOS2_w);
    }

    for (j = 0; j < 3; j++)
      force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
    if (fac * dist > 1000)
      fprintf(stderr, "%d: LJ-Warning: Pair (%d-%d) force=%f dist=%f\n",
              this_node, p1->p.identity, p2->p.identity, fac * dist, dist);
#endif

    ONEPART_TRACE(if (p1->p.identity == check_id)
                      fprintf(stderr, "%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with "
                                      "part id=%d at dist %f fac %.3e\n",
                              this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                              p2->p.identity, dist, fac));
    ONEPART_TRACE(if (p2->p.identity == check_id)
                      fprintf(stderr, "%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with "
                                      "part id=%d at dist %f fac %.3e\n",
                              this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                              p1->p.identity, dist, fac));

    LJ_TRACE(fprintf(
        stderr, "%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
        this_node, p1->p.identity, p2->p.identity, dist, fac * d[0], fac * d[1],
        fac * d[2]));
  }
}

/** calculate lj-cos2 energy between particle p1 and p2. */
inline double ljcos2_pair_energy(Particle *p1, Particle *p2,
                                 IA_parameters *ia_params, double d[3],
                                 double dist) {
  double r_off, frac2, frac6;

  if ((dist < ia_params->LJCOS2_cut + ia_params->LJCOS2_offset)) {
    r_off = dist - ia_params->LJCOS2_offset;
    if (r_off < ia_params->LJCOS2_rchange) {
      frac2 = SQR(ia_params->LJCOS2_sig / r_off);
      frac6 = frac2 * frac2 * frac2;
      return 4.0 * ia_params->LJCOS2_eps * (SQR(frac6) - frac6);
    } else if (r_off < ia_params->LJCOS2_rchange + ia_params->LJCOS2_w) {
      return -ia_params->LJCOS2_eps / 2 *
             (cos(M_PI * (r_off - ia_params->LJCOS2_rchange) /
                  ia_params->LJCOS2_w) +
              1);
    }
  }
  return 0.0;
}

#endif /* ifdef LJCOS2 */
#endif
