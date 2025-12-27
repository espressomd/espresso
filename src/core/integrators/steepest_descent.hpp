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

#pragma once

#include "cell_system/CellStructure.hpp"

/** @brief Steepest descent algorithm. */
struct SteepestDescent {
  /**
   * @brief Maximal particle force or torque.
   *
   * If the maximal force experienced by particles in the system (in any
   * direction) is inferior to this threshold, minimization stops.
   */
  double f_max = 0.;
  /** Dampening constant */
  double gamma = 0.;
  /**
   * Maximal particle displacement or rotation.
   *
   * Maximal distance in MD units of length or rotation in radians that
   * a particle can experience during one integration step, in one direction.
   */
  double max_displacement = 0.;

  SteepestDescent() = default;
  SteepestDescent(double f_max, double gamma, double max_displacement);

  /**
   * @brief Run steepest descent algorithm.
   * @return whether the maximum force/torque encountered is below the user
   *         limit @ref SteepestDescent::f_max "f_max".
   */
  bool propagate(CellStructure &cell_structure) const;
};
