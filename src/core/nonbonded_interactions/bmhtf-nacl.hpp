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
#ifndef BMHTF_NACL_H
#define BMHTF_NACL_H
/** \file
 *  Routines to calculate the Born-Meyer-Huggins-Tosi-Fumi potential
 *  between particle pairs.
 *
 *  Implementation in \ref bmhtf-nacl.cpp.
 */

#include "nonbonded_interaction_data.hpp"

#include <utils/math/int_pow.hpp>

#ifdef BMHTF_NACL

int BMHTF_set_params(int part_type_a, int part_type_b, double A, double B,
                     double C, double D, double sig, double cut);

/** Calculate BMHTF force factor */
inline double BMHTF_pair_force_factor(IA_parameters const &ia_params,
                                      double dist) {
  if (dist < ia_params.bmhtf.cut) {
    auto const dist8 = Utils::int_pow<8>(dist);
    auto const dist10 = Utils::int_pow<10>(dist);
    return ia_params.bmhtf.A * ia_params.bmhtf.B *
               exp(ia_params.bmhtf.B * (ia_params.bmhtf.sig - dist)) / dist -
           6 * ia_params.bmhtf.C / dist8 - 8 * ia_params.bmhtf.D / dist10;
  }
  return 0.0;
}

/** Calculate BMHTF force */
inline Utils::Vector3d BMHTF_pair_force(IA_parameters const &ia_params,
                                        Utils::Vector3d const &d, double dist) {
  return d * BMHTF_pair_force_factor(ia_params, dist);
}

/** Calculate BMHTF potential energy */
inline double BMHTF_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < ia_params.bmhtf.cut) {
    auto const dist6 = Utils::int_pow<6>(dist);
    auto const dist8 = Utils::int_pow<8>(dist);
    return ia_params.bmhtf.A *
               exp(ia_params.bmhtf.B * (ia_params.bmhtf.sig - dist)) -
           ia_params.bmhtf.C / dist6 - ia_params.bmhtf.D / dist8 +
           ia_params.bmhtf.computed_shift;
  }
  return 0.0;
}

#endif
#endif
