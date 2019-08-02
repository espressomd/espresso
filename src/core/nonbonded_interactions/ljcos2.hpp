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
#ifndef _LJCOS2_H
#define _LJCOS2_H

/** \file
 *  Routines to calculate the Lennard-Jones with cosine tail potential
 *  between particle pairs. Cosine tail is different from that in
 *  \ref ljcos.hpp. Used for attractive tail/tail interactions in lipid
 *  bilayer calculations.
 *
 *  Implementation in \ref ljcos2.cpp.
 */

#include "config.hpp"

#ifdef LJCOS2

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <cmath>
#include <utils/math/int_pow.hpp>

int ljcos2_set_params(int part_type_a, int part_type_b, double eps, double sig,
                      double offset, double w);

/** Calculate lj-cos2 force between particle p1 and p2. */
inline void add_ljcos2_pair_force(Particle const *const p1,
                                  Particle const *const p2,
                                  IA_parameters const *const ia_params,
                                  Utils::Vector3d const &d, double dist,
                                  Utils::Vector3d &force) {
  if (dist < (ia_params->ljcos2.cut + ia_params->ljcos2.offset)) {
    auto const r_off = dist - ia_params->ljcos2.offset;
    auto fac = 0.0;
    if (r_off < ia_params->ljcos2.rchange) {
      auto const frac6 = Utils::int_pow<6>(ia_params->ljcos2.sig / r_off);
      fac =
          48.0 * ia_params->ljcos2.eps * frac6 * (frac6 - 0.5) / (r_off * dist);
    } else if (r_off < ia_params->ljcos2.rchange + ia_params->ljcos2.w) {
      fac =
          -ia_params->ljcos2.eps * M_PI / 2 / ia_params->ljcos2.w / dist *
          sin(M_PI * (r_off - ia_params->ljcos2.rchange) / ia_params->ljcos2.w);
    }
    force += fac * d;

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
        this_node, p1->p.identity, p2->p.identity, dist, force[0], force[1],
        force[2]));
  }
}

/** Calculate lj-cos2 energy between particle p1 and p2. */
inline double ljcos2_pair_energy(Particle const *const p1,
                                 Particle const *const p2,
                                 IA_parameters const *const ia_params,
                                 Utils::Vector3d const &d, double dist) {
  if (dist < (ia_params->ljcos2.cut + ia_params->ljcos2.offset)) {
    auto const r_off = dist - ia_params->ljcos2.offset;
    if (r_off < ia_params->ljcos2.rchange) {
      auto const frac6 = Utils::int_pow<6>(ia_params->ljcos2.sig / r_off);
      return 4.0 * ia_params->ljcos2.eps * (Utils::sqr(frac6) - frac6);
    }
    if (r_off < (ia_params->ljcos2.rchange + ia_params->ljcos2.w)) {
      auto const fac = -ia_params->ljcos2.eps / 2 *
                       (cos(M_PI * (r_off - ia_params->ljcos2.rchange) /
                            ia_params->ljcos2.w) +
                        1);
      return fac;
    }
  }
  return 0.0;
}

#endif /* ifdef LJCOS2 */
#endif
