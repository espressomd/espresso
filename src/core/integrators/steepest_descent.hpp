/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef CORE_INTEGRATORS_STEEPEST_DESCENT_HPP
#define CORE_INTEGRATORS_STEEPEST_DESCENT_HPP

#include "ParticleRange.hpp"

/** Parameters for the steepest descent algorithm */
struct SteepestDescentParameters {
  /** Maximal particle force
   *
   *  If the maximal force experienced by particles in the system (in any
   *  direction) is inferior to this threshold, minimization stops.
   */
  double f_max;
  /** Dampening constant */
  double gamma;
  /** Maximal particle displacement
   *
   *  Maximal distance that a particle can travel during one integration step,
   *  in one direction.
   */
  double max_displacement;

  SteepestDescentParameters(double f_max, double gamma,
                            double max_displacement);
};

void register_integrator(SteepestDescentParameters const &obj);

/** Steepest descent integrator
 *  @return whether the maximum force/torque encountered is below the user
 *          limit @ref SteepestDescentParameters::f_max "f_max".
 */
bool steepest_descent_step(const ParticleRange &particles);

#endif /* CORE_INTEGRATORS_STEEPEST_DESCENT_HPP */
