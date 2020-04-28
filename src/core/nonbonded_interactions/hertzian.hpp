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
#ifndef HERTZIAN_H
#define HERTZIAN_H

/** \file
 *  Routines to calculate the Hertzian potential between particle pairs.
 *
 *  Implementation in \ref hertzian.cpp.
 */

#include "nonbonded_interaction_data.hpp"

#ifdef HERTZIAN

int hertzian_set_params(int part_type_a, int part_type_b, double eps,
                        double sig);

/** Calculate Hertzian force factor */
inline double hertzian_pair_force_factor(IA_parameters const &ia_params,
                                         double dist) {
  if (dist < ia_params.hertzian.sig) {
    return 5. / 2. * ia_params.hertzian.eps / ia_params.hertzian.sig *
           pow(1 - dist / ia_params.hertzian.sig, 3. / 2.) / dist;
  }
  return 0.0;
}

/** Calculate Hertzian force */
inline Utils::Vector3d hertzian_pair_force(IA_parameters const &ia_params,
                                           Utils::Vector3d const &d,
                                           double dist) {
  return d * hertzian_pair_force_factor(ia_params, dist);
}

/** Calculate Hertzian energy */
inline double hertzian_pair_energy(IA_parameters const &ia_params,
                                   double dist) {
  if (dist < ia_params.hertzian.sig) {
    return ia_params.hertzian.eps *
           pow(1 - dist / ia_params.hertzian.sig, 5. / 2.);
  }
  return 0.0;
}

#endif /* ifdef HERTZIAN */
#endif
