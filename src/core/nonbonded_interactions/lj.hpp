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
#ifndef _LJ_H
#define _LJ_H

#include "config.hpp"

#ifdef LENNARD_JONES

/** \file
 *  Routines to calculate the Lennard-Jones potential between particle pairs.
 *
 *  Implementation in \ref lj.cpp.
 */

#include "nonbonded_interaction_data.hpp"

#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

int lennard_jones_set_params(int part_type_a, int part_type_b, double eps,
                             double sig, double cut, double shift,
                             double offset, double min);

/** Calculate Lennard-Jones force factor */
inline double lj_pair_force_factor(IA_parameters const &ia_params,
                                   double dist) {
  if ((dist < ia_params.lj.cut + ia_params.lj.offset) &&
      (dist > ia_params.lj.min + ia_params.lj.offset)) {
    auto const r_off = dist - ia_params.lj.offset;
    auto const frac6 = Utils::int_pow<6>(ia_params.lj.sig / r_off);
    return 48.0 * ia_params.lj.eps * frac6 * (frac6 - 0.5) / (r_off * dist);
  }
  return 0.0;
}

/** Calculate Lennard-Jones force */
inline Utils::Vector3d lj_pair_force(IA_parameters const &ia_params,
                                     Utils::Vector3d const &d, double dist) {
  return d * lj_pair_force_factor(ia_params, dist);
}

/** Calculate Lennard-Jones energy */
inline double lj_pair_energy(IA_parameters const &ia_params, double dist) {
  if ((dist < ia_params.lj.cut + ia_params.lj.offset) &&
      (dist > ia_params.lj.min + ia_params.lj.offset)) {
    auto const r_off = dist - ia_params.lj.offset;
    auto const frac6 = Utils::int_pow<6>(ia_params.lj.sig / r_off);
    return 4.0 * ia_params.lj.eps *
           (Utils::sqr(frac6) - frac6 + ia_params.lj.shift);
  }
  return 0.0;
}

#endif /* ifdef LENNARD_JONES */
#endif
