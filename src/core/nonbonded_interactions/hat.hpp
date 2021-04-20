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
#ifndef hat_H
#define hat_H

/** \file
 *  Routines to calculate the hat potential between particle pairs.
 *
 *  Implementation in \ref hat.cpp.
 */

#include "config.hpp"

#ifdef HAT

#include "nonbonded_interaction_data.hpp"

#include <utils/Vector.hpp>

int hat_set_params(int part_type_a, int part_type_b, double Fmax, double r);

/** Calculate hat force factor */
inline double hat_pair_force_factor(IA_parameters const &ia_params,
                                    double dist) {
  if (dist > 0. && dist < ia_params.hat.r) {
    return ia_params.hat.Fmax * (1.0 - dist / ia_params.hat.r) / dist;
  }
  return 0.0;
}

/** Calculate hat energy */
inline double hat_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < ia_params.hat.r) {
    auto const r = ia_params.hat.r;
    return ia_params.hat.Fmax * (dist - r) * ((dist + r) / (2.0 * r) - 1.0);
  }
  return 0.0;
}

#endif /* ifdef HAT */
#endif
