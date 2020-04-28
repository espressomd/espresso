/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef CORE_TABULATED_POTENTIAL_HPP
#define CORE_TABULATED_POTENTIAL_HPP

#include <utils/linear_interpolation.hpp>

#include <boost/algorithm/clamp.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <vector>

/** Evaluate forces and energies using a custom potential profile.
 *
 *  Forces and energies are evaluated by linear interpolation. The curves
 *  @ref force_tab and @ref energy_tab must be sampled uniformly between
 *  @ref minval and @ref maxval.
 */
struct TabulatedPotential {
  /** Position on the x-axis of the first tabulated value. */
  double minval = -1.0;
  /** Position on the x-axis of the last tabulated value. */
  double maxval = -1.0;
  /** %Distance on the x-axis between tabulated values. */
  double invstepsize = 0.0;
  /** Tabulated forces. */
  std::vector<double> force_tab;
  /** Tabulated energies. */
  std::vector<double> energy_tab;

  /** Evaluate the force at position @p x.
   *  @param x  Bond length/angle
   *  @return Interpolated force.
   */
  double force(double x) const {
    using boost::algorithm::clamp;
    return Utils::linear_interpolation(force_tab, invstepsize, minval,
                                       clamp(x, minval, maxval));
  }

  /** Evaluate the energy at position @p x.
   *  @param x  Bond length/angle
   *  @return Interpolated energy.
   */
  double energy(double x) const {
    using boost::algorithm::clamp;
    return Utils::linear_interpolation(energy_tab, invstepsize, minval,
                                       clamp(x, minval, maxval));
  }

  double cutoff() const { return maxval; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &minval;
    ar &maxval;
    ar &invstepsize;
    ar &force_tab;
    ar &energy_tab;
  }
};

#endif
