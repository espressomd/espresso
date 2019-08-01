/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef STEPPOT_H
#define STEPPOT_H

/** \file
 *  Routines to calculate the smooth step potential between particle pairs.
 *
 *  Implementation in \ref steppot.cpp.
 */

#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#ifdef SMOOTH_STEP

int smooth_step_set_params(int part_type_a, int part_type_b, double d, int n,
                           double eps, double k0, double sig, double cut);

/** Calculate smooth step force between particle p1 and p2 */
inline Utils::Vector3d add_SmSt_pair_force(Particle const *const p1,
                                           Particle const *const p2,
                                           IA_parameters const *const ia_params,
                                           Utils::Vector3d const &d,
                                           double dist, double dist2) {
  if (dist >= ia_params->SmSt_cut) {
    return {};
  }

  auto const frac = ia_params->SmSt_d / dist;
  auto const fracP = pow(frac, ia_params->SmSt_n);
  auto const er = exp(2. * ia_params->SmSt_k0 * (dist - ia_params->SmSt_sig));
  auto const fac = (ia_params->SmSt_n * fracP + 2. * ia_params->SmSt_eps *
                                                    ia_params->SmSt_k0 * dist *
                                                    er / Utils::sqr(1.0 + er)) /
                   dist2;
  auto const force = fac * d;
  return force;
}

/** Calculate smooth step energy between particle p1 and p2. */
inline double SmSt_pair_energy(Particle const *const p1,
                               Particle const *const p2,
                               IA_parameters const *const ia_params,
                               Utils::Vector3d const &d, double dist,
                               double dist2) {
  if (dist >= ia_params->SmSt_cut) {
    return 0.0;
  }

  auto const frac = ia_params->SmSt_d / dist;
  auto const fracP = pow(frac, ia_params->SmSt_n);
  auto const er = exp(2. * ia_params->SmSt_k0 * (dist - ia_params->SmSt_sig));
  auto const fac = fracP + ia_params->SmSt_eps / (1.0 + er);

  return fac;
}

#endif /* ifdef SMOOTH_STEP */
#endif
