/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _LJGEN_H
#define _LJGEN_H

#include "config.hpp"

#ifdef LENNARD_JONES_GENERIC

/** \file
 *  Routines to calculate the generalized Lennard-Jones potential between
 *  particle pairs. "Generalized" here means that the LJ energy is of the form
 *  @f[
 *      \varepsilon \cdot
 *      \left[
 *            b_1 \left(\frac{\sigma}{r-r_{\text{offset}}}\right)^{a_1}
 *          - b_2 \left(\frac{\sigma}{r-r_{\text{offset}}}\right)^{a_2}
 *          + \text{shift}
 *      \right]
 *  @f]
 *
 *  Implementation in \ref ljgen.cpp.
 */

#include "nonbonded_interaction_data.hpp"

int ljgen_set_params(int part_type_a, int part_type_b, double eps, double sig,
                     double cut, double shift, double offset, double a1,
                     double a2, double b1, double b2
#ifdef LJGEN_SOFTCORE
                     ,
                     double lambda, double softrad
#endif
);

/** Calculate Lennard-Jones force factor */
inline double ljgen_pair_force_factor(IA_parameters const &ia_params,
                                      double dist) {
  if (dist < (ia_params.ljgen.cut + ia_params.ljgen.offset)) {
    auto r_off = dist - ia_params.ljgen.offset;

#ifdef LJGEN_SOFTCORE
    r_off *= r_off;
    r_off += Utils::sqr(ia_params.ljgen.sig) * (1.0 - ia_params.ljgen.lambda1) *
             ia_params.ljgen.softrad;
    r_off = sqrt(r_off);
#else
    r_off = std::abs(r_off);
#endif
    auto const frac = ia_params.ljgen.sig / r_off;
    auto const fac = ia_params.ljgen.eps
#ifdef LJGEN_SOFTCORE
                     * ia_params.ljgen.lambda1 *
                     (dist - ia_params.ljgen.offset) / r_off
#endif
                     * (ia_params.ljgen.b1 * ia_params.ljgen.a1 *
                            pow(frac, ia_params.ljgen.a1) -
                        ia_params.ljgen.b2 * ia_params.ljgen.a2 *
                            pow(frac, ia_params.ljgen.a2)) /
                     (r_off * dist);
    return fac;
  }
  return 0.0;
}

/** Calculate Lennard-Jones force */
inline Utils::Vector3d ljgen_pair_force(IA_parameters const &ia_params,
                                        Utils::Vector3d const &d, double dist) {
  return d * ljgen_pair_force_factor(ia_params, dist);
}

/** Calculate Lennard-Jones energy */
inline double ljgen_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < (ia_params.ljgen.cut + ia_params.ljgen.offset)) {
    auto r_off = dist - ia_params.ljgen.offset;
#ifdef LJGEN_SOFTCORE
    r_off *= r_off;
    r_off += pow(ia_params.ljgen.sig, 2) * (1.0 - ia_params.ljgen.lambda1) *
             ia_params.ljgen.softrad;
    r_off = sqrt(r_off);
#else
    r_off = std::abs(r_off);
#endif
    auto const frac = ia_params.ljgen.sig / r_off;
    auto const fac = ia_params.ljgen.eps
#ifdef LJGEN_SOFTCORE
                     * ia_params.ljgen.lambda1
#endif
                     * (ia_params.ljgen.b1 * pow(frac, ia_params.ljgen.a1) -
                        ia_params.ljgen.b2 * pow(frac, ia_params.ljgen.a2) +
                        ia_params.ljgen.shift);
    return fac;
  }
  return 0.0;
}

#endif

/* LJGEN_H */
#endif
