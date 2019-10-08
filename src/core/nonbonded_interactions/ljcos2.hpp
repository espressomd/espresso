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

#include "nonbonded_interaction_data.hpp"

#include <cmath>
#include <utils/math/int_pow.hpp>

int ljcos2_set_params(int part_type_a, int part_type_b, double eps, double sig,
                      double offset, double w);

/** Calculate Lennard-Jones cosine squared force factor */
inline double ljcos2_pair_force_factor(IA_parameters const &ia_params,
                                       double dist) {
  if (dist < (ia_params.ljcos2.cut + ia_params.ljcos2.offset)) {
    auto const r_off = dist - ia_params.ljcos2.offset;
    auto fac = 0.0;
    if (r_off < ia_params.ljcos2.rchange) {
      auto const frac6 = Utils::int_pow<6>(ia_params.ljcos2.sig / r_off);
      fac =
          48.0 * ia_params.ljcos2.eps * frac6 * (frac6 - 0.5) / (r_off * dist);
    } else if (r_off < ia_params.ljcos2.rchange + ia_params.ljcos2.w) {
      fac = -ia_params.ljcos2.eps * M_PI / 2 / ia_params.ljcos2.w / dist *
            sin(M_PI * (r_off - ia_params.ljcos2.rchange) / ia_params.ljcos2.w);
    }
    return fac;
  }
  return 0.0;
}

/** Calculate Lennard-Jones cosine squared force */
inline Utils::Vector3d ljcos2_pair_force(IA_parameters const &ia_params,
                                         Utils::Vector3d const &d,
                                         double dist) {
  return d * ljcos2_pair_force_factor(ia_params, dist);
}

/** Calculate Lennard-Jones cosine squared energy */
inline double ljcos2_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < (ia_params.ljcos2.cut + ia_params.ljcos2.offset)) {
    auto const r_off = dist - ia_params.ljcos2.offset;
    if (r_off < ia_params.ljcos2.rchange) {
      auto const frac6 = Utils::int_pow<6>(ia_params.ljcos2.sig / r_off);
      return 4.0 * ia_params.ljcos2.eps * (Utils::sqr(frac6) - frac6);
    }
    if (r_off < (ia_params.ljcos2.rchange + ia_params.ljcos2.w)) {
      return -ia_params.ljcos2.eps / 2 *
             (cos(M_PI * (r_off - ia_params.ljcos2.rchange) /
                  ia_params.ljcos2.w) +
              1);
    }
  }
  return 0.0;
}

#endif /* ifdef LJCOS2 */
#endif
