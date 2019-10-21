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

#ifndef __MINIMIZE_ENERGY_HPP
#define __MINIMIZE_ENERGY_HPP

#include "ParticleRange.hpp"

/** Parameters for the steepest descent algorithm */
struct MinimizeEnergyParameters {
  /** Maximal particle force
   *
   *  If the maximal force experienced by particles in the system (in any
   *  direction) is inferior to this threshold, minimization stops.
   */
  double f_max;
  /** Dampening constant */
  double gamma;
  /** Maximal number of integration steps */
  int max_steps;
  /** Maximal particle displacement
   *
   *  Maximal distance that a particle can travel during one integration step,
   *  in one direction.
   */
  double max_displacement;
};

/** Steepest descent initializer
 *
 *  Sets the parameters in @ref MinimizeEnergyParameters
 */
void minimize_energy_init(double f_max, double gamma, int max_steps,
                          double max_displacement);

/** Steepest descent main integration loop
 *
 *  Integration stops when the maximal force is lower than the user limit
 *  @ref MinimizeEnergyParameters::f_max "f_max" or when the maximal number
 *  of steps @ref MinimizeEnergyParameters::max_steps "max_steps" is reached.
 */
void minimize_energy();

/** Steepest descent integrator
 *  @return whether the maximum force/torque encountered is below the user
 *          limit @ref MinimizeEnergyParameters::f_max "f_max".
 */
bool steepest_descent_step(const ParticleRange &particles);

#endif /* __MINIMIZE_ENERGY */
