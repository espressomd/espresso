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
#ifndef DEBYE_HUECKEL_H
#define DEBYE_HUECKEL_H
/** \file
 *  Routines to calculate the Debye-Hückel energy and force
 *  for a particle pair.
 */
#include "config.hpp"

#ifdef ELECTROSTATICS

#include <utils/Vector.hpp>

#include <cmath>

/** Debye-Hückel parameters. */
struct Debye_hueckel_params {
  /** Interaction cutoff. */
  double r_cut;
  /** Ionic strength. */
  double kappa;
};

/** Global state of the Debye-Hückel method. */
extern Debye_hueckel_params dh_params;

void dh_set_params(double kappa, double r_cut);

/** Compute the Debye-Hueckel pair force.
 *  @param[in]  q1q2      Product of the charges on p1 and p2.
 *  @param[in]  d         Vector pointing from p1 to p2.
 *  @param[in]  dist      Distance between p1 and p2.
 *  @param[out] force     Calculated force on p1.
 */
inline void add_dh_coulomb_pair_force(double const q1q2,
                                      Utils::Vector3d const &d,
                                      double const dist,
                                      Utils::Vector3d &force) {
  if (dist < dh_params.r_cut) {
    // pure Coulomb case
    auto fac = q1q2 / (dist * dist * dist);
    if (dh_params.kappa > 0.0) {
      // Debye-Hueckel case
      auto const kappa_dist = dh_params.kappa * dist;
      fac *= exp(-kappa_dist) * (1.0 + kappa_dist);
    }
    force += fac * d;
  }
}

/** Compute the Debye-Hueckel pair energy.
 *  @param q1q2      Product of the charges on p1 and p2.
 *  @param dist      Distance between p1 and p2.
 */
inline double dh_coulomb_pair_energy(double const q1q2, double const dist) {
  if (dist < dh_params.r_cut) {
    if (dh_params.kappa > 0.0)
      return q1q2 * exp(-dh_params.kappa * dist) / dist;

    return q1q2 / dist;
  }
  return 0.0;
}
#endif

#endif
