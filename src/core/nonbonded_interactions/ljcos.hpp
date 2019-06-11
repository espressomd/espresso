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
#ifndef _LJCOS_H
#define _LJCOS_H
/** \file
 *  Routines to calculate the Lennard-Jones+cosine energy and/or force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "config.hpp"

#ifdef LJCOS

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

int ljcos_set_params(int part_type_a, int part_type_b, double eps, double sig,
                     double cut, double offset);

inline void add_ljcos_pair_force(const Particle *const p1,
                                 const Particle *const p2,
                                 IA_parameters *ia_params, double const d[3],
                                 double dist, double force[3]) {
  if ((dist < ia_params->LJCOS_cut + ia_params->LJCOS_offset)) {
    double r_off = dist - ia_params->LJCOS_offset;
    /* cos part of ljcos potential. */
    if (dist > ia_params->LJCOS_rmin + ia_params->LJCOS_offset) {
      double fac = (r_off / dist) * ia_params->LJCOS_alfa *
                   ia_params->LJCOS_eps *
                   (sin(ia_params->LJCOS_alfa * Utils::sqr(r_off) +
                        ia_params->LJCOS_beta));
      for (int j = 0; j < 3; j++)
        force[j] += fac * d[j];
    }
    /* Lennard-Jones part of the potential. */
    else if (dist > 0) {
      double frac2 = Utils::sqr(ia_params->LJCOS_sig / r_off);
      double frac6 = frac2 * frac2 * frac2;
      double fac =
          48.0 * ia_params->LJCOS_eps * frac6 * (frac6 - 0.5) / (r_off * dist);

      for (int j = 0; j < 3; j++)
        force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if (fac * dist > 1000)
        fprintf(stderr, "%d: LJCOS-Warning: Pair (%d-%d) force=%f dist=%f\n",
                this_node, p1->p.identity, p2->p.identity, fac * dist, dist);
#endif
    }
  }
}

inline double ljcos_pair_energy(const Particle *p1, const Particle *p2,
                                const IA_parameters *ia_params,
                                const double d[3], double dist) {
  if ((dist < ia_params->LJCOS_cut + ia_params->LJCOS_offset)) {
    double r_off = dist - ia_params->LJCOS_offset;
    /* Lennard-Jones part of the potential. */
    if (dist < (ia_params->LJCOS_rmin + ia_params->LJCOS_offset)) {
      double frac2 = Utils::sqr(ia_params->LJCOS_sig / r_off);
      double frac6 = frac2 * frac2 * frac2;
      return 4.0 * ia_params->LJCOS_eps * (Utils::sqr(frac6) - frac6);
    }
    /* cosine part of the potential. */
    if (dist < (ia_params->LJCOS_cut + ia_params->LJCOS_offset)) {
      return .5 * ia_params->LJCOS_eps *
             (cos(ia_params->LJCOS_alfa * Utils::sqr(r_off) +
                  ia_params->LJCOS_beta) -
              1.);
    }
    /* this should not happen! */

    fprintf(stderr, "this is the distance, which is negative %.3e\n", r_off);
  }
  return 0.0;
}

#endif
#endif
