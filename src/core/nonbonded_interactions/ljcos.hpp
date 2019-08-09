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
 *  Routines to calculate the Lennard-Jones+cosine potential between
 *  particle pairs.
 *
 *  Implementation in \ref ljcos.cpp.
 */

#include "config.hpp"

#ifdef LJCOS

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <utils/math/int_pow.hpp>

int ljcos_set_params(int part_type_a, int part_type_b, double eps, double sig,
                     double cut, double offset);

inline void add_ljcos_pair_force(Particle const *const p1,
                                 Particle const *const p2,
                                 IA_parameters const *const ia_params,
                                 Utils::Vector3d const &d, double dist,
                                 Utils::Vector3d &force) {
  if ((dist < ia_params->ljcos.cut + ia_params->ljcos.offset)) {
    auto const r_off = dist - ia_params->ljcos.offset;
    /* cos part of ljcos potential. */
    if (dist > ia_params->ljcos.rmin + ia_params->ljcos.offset) {
      auto const fac = (r_off / dist) * ia_params->ljcos.alfa *
                       ia_params->ljcos.eps *
                       (sin(ia_params->ljcos.alfa * Utils::sqr(r_off) +
                            ia_params->ljcos.beta));
      force += fac * d;
    }
    /* Lennard-Jones part of the potential. */
    else if (dist > 0) {
      auto const frac6 = Utils::int_pow<6>(ia_params->ljcos.sig / r_off);
      auto const fac =
          48.0 * ia_params->ljcos.eps * frac6 * (frac6 - 0.5) / (r_off * dist);
      force += fac * d;

#ifdef LJ_WARN_WHEN_CLOSE
      if (fac * dist > 1000)
        fprintf(stderr, "%d: LJCOS-Warning: Pair (%d-%d) force=%f dist=%f\n",
                this_node, p1->p.identity, p2->p.identity, fac * dist, dist);
#endif
    }
  }
}

inline double ljcos_pair_energy(Particle const *const p1,
                                Particle const *const p2,
                                IA_parameters const *const ia_params,
                                Utils::Vector3d const &d, double dist) {
  if (dist < (ia_params->ljcos.cut + ia_params->ljcos.offset)) {
    auto const r_off = dist - ia_params->ljcos.offset;
    /* Lennard-Jones part of the potential. */
    if (dist < (ia_params->ljcos.rmin + ia_params->ljcos.offset)) {
      auto const frac6 = Utils::int_pow<6>(ia_params->ljcos.sig / r_off);
      return 4.0 * ia_params->ljcos.eps * (Utils::sqr(frac6) - frac6);
    }
    /* cosine part of the potential. */
    if (dist < (ia_params->ljcos.cut + ia_params->ljcos.offset)) {
      auto const fac = .5 * ia_params->ljcos.eps *
                       (cos(ia_params->ljcos.alfa * Utils::sqr(r_off) +
                            ia_params->ljcos.beta) -
                        1.);
      return fac;
    }
    /* this should not happen! */
    fprintf(stderr, "this is the distance, which is negative %.3e\n", r_off);
  }
  return 0.0;
}

#endif
#endif
