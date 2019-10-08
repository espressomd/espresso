/*
 * Copyright (C) 2018-2019 The ESPResSo project
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
#ifndef WCA_HPP
#define WCA_HPP
/** \file
 *  Routines to calculate the Weeks-Chandler-Andersen potential between
 *  particle pairs.
 *
 *  Implementation in \ref wca.cpp.
 */

#include "config.hpp"

#ifdef WCA

#include "nonbonded_interaction_data.hpp"
#include <utils/math/int_pow.hpp>

int wca_set_params(int part_type_a, int part_type_b, double eps, double sig);

/** Calculate WCA force factor */
inline double wca_pair_force_factor(IA_parameters const &ia_params,
                                    double dist) {
  if (dist < ia_params.wca.cut) {
    auto const frac6 = Utils::int_pow<6>(ia_params.wca.sig / dist);
    return 48.0 * ia_params.wca.eps * frac6 * (frac6 - 0.5) / (dist * dist);
  }
  return 0.0;
}

/** Calculate WCA force */
inline Utils::Vector3d wca_pair_force(IA_parameters const &ia_params,
                                      Utils::Vector3d const &d, double dist) {
  return d * wca_pair_force_factor(ia_params, dist);
}

/** Calculate WCA energy */
inline double wca_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < ia_params.wca.cut) {
    auto const frac6 = Utils::int_pow<6>(ia_params.wca.sig / dist);
    return 4.0 * ia_params.wca.eps * (Utils::sqr(frac6) - frac6 + .25);
  }
  return 0.0;
}

#endif /* ifdef WCA */
#endif
