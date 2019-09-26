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
#ifndef soft_H
#define soft_H

/** \file
 *  Routines to calculate the soft-sphere potential between particle pairs.
 *
 *  Implementation in \ref soft_sphere.cpp
 */

#include "config.hpp"

#ifdef SOFT_SPHERE

#include "nonbonded_interaction_data.hpp"

int soft_sphere_set_params(int part_type_a, int part_type_b, double a, double n,
                           double cut, double offset);

/** Resultant force due to a soft-sphere potential between two
 *  particles at interatomic separation r
 */
inline double soft_force_r(double a, double n, double r) {
  return (a * n / pow(r, n + 1));
}

/** Potential energy due to a soft-sphere potential between two
 *  particles at interatomic separation r
 */
inline double soft_energy_r(double a, double n, double r) {
  return (a / pow(r, n));
}

/** Calculate soft-sphere force factor */
inline double soft_pair_force_factor(IA_parameters const &ia_params,
                                     double dist) {
  if (dist < (ia_params.soft_sphere.cut + ia_params.soft_sphere.offset)) {
    /* normal case: resulting force/energy smaller than zero. */
    auto const r_off = dist - ia_params.soft_sphere.offset;
    if (r_off > 0.0) {
      return soft_force_r(ia_params.soft_sphere.a, ia_params.soft_sphere.n,
                          r_off) /
             dist;
    }
  }
  return 0.0;
}

/** Calculate soft-sphere force */
inline Utils::Vector3d soft_pair_force(IA_parameters const &ia_params,
                                       Utils::Vector3d const &d, double dist) {
  return d * soft_pair_force_factor(ia_params, dist);
}

/** Calculate soft-sphere energy */
inline double soft_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < (ia_params.soft_sphere.cut + ia_params.soft_sphere.offset)) {
    auto const r_off = dist - ia_params.soft_sphere.offset;
    /* normal case: resulting force/energy smaller than zero. */
    return soft_energy_r(ia_params.soft_sphere.a, ia_params.soft_sphere.n,
                         r_off);
  }
  return 0.0;
}

#endif /* ifdef SOFT_SPHERE */
#endif
