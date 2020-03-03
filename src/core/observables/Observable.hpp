/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#ifndef OBSERVABLES_OBSERVABLE_HPP
#define OBSERVABLES_OBSERVABLE_HPP

#include <fstream>
#include <numeric>
#include <string>
#include <vector>

#include "PartCfg.hpp"

namespace Observables {

/** Base class for observables.
 *
 *  An observable extracts raw data from a system or compute a statistic based
 *  on the state of a system, and returns an array of doubles.
 *
 *  %Observables typically don't have setters or getters to access and modify
 *  their member variables, and usually only have a default constructor with no
 *  argument. Each observable class has a corresponding interface in
 *  @ref ScriptInterface::Observables, where setters and getters are defined.
 */
class Observable {
public:
  Observable() = default;
  virtual ~Observable() = default;
  /** Calculate the set of values measured by the observable */
  virtual std::vector<double> operator()() const = 0;

  /** Size of the flat array returned by the observable */
  int n_values() const {
    auto const v = shape();
    return static_cast<int>(
        std::accumulate(v.begin(), v.end(), 1, std::multiplies<>()));
  }

  /** Dimensions needed to reshape the flat array returned by the observable */
  virtual std::vector<size_t> shape() const = 0;
};

} // Namespace Observables
#endif
