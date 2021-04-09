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
#ifndef REACTION_FIELD_H
#define REACTION_FIELD_H
/** \file
 *  Routines to calculate the Reaction Field energy or/and force
 *  for a particle pair @cite neumann85b, @cite tironi95a.
 *
 *  Implementation in \ref reaction_field.cpp
 */

#include "config.hpp"

#ifdef ELECTROSTATICS

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

/** Reaction Field parameters. */
struct Reaction_field_params {
  /** Ionic strength. */
  double kappa;
  /** Continuum dielectric constant inside the cavity. */
  double epsilon1;
  /** Continuum dielectric constant outside the cavity. */
  double epsilon2;
  /** Interaction cutoff. */
  double r_cut;
  /** Interaction prefactor. Corresponds to the quantity
   *  @f$ 1 + B_1 @f$ from eq. 22 in @cite tironi95a.
   */
  double B;
};

/** Global state of the Reaction Field method. */
extern Reaction_field_params rf_params;

void rf_set_params(double kappa, double epsilon1, double epsilon2,
                   double r_cut);

/** Compute the Reaction Field pair force.
 *  @param[in]  q1q2      Product of the charges on p1 and p2.
 *  @param[in]  d         Vector pointing from p1 to p2.
 *  @param[in]  dist      Distance between p1 and p2.
 *  @param[out] force     Calculated force on p1.
 */
inline void add_rf_coulomb_pair_force(double const q1q2,
                                      Utils::Vector3d const &d,
                                      double const dist,
                                      Utils::Vector3d &force) {
  if (dist < rf_params.r_cut) {
    auto fac = 1.0 / Utils::int_pow<3>(dist) +
               rf_params.B / Utils::int_pow<3>(rf_params.r_cut);
    fac *= q1q2;
    force += fac * d;
  }
}

/** Compute the Reaction Field pair energy.
 *  @param q1q2      Product of the charges on p1 and p2.
 *  @param dist      Distance between p1 and p2.
 */
inline double rf_coulomb_pair_energy(double const q1q2, double const dist) {
  if (dist < rf_params.r_cut) {
    auto fac = 1.0 / dist - (rf_params.B * dist * dist) /
                                (2 * Utils::int_pow<3>(rf_params.r_cut));
    // remove discontinuity at dist = r_cut
    fac -= (1 - rf_params.B / 2) / rf_params.r_cut;
    fac *= q1q2;
    return fac;
  }
  return 0.0;
}

#endif

#endif
