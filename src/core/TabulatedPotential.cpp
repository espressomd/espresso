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

#include "TabulatedPotential.hpp"

#include <stdexcept>
#include <vector>

TabulatedPotential::TabulatedPotential(double minval, double maxval,
                                       std::vector<double> const &force,
                                       std::vector<double> const &energy)
    : minval{minval}, maxval{maxval} {

  if (minval > maxval) {
    throw std::domain_error("TabulatedPotential parameter 'max' must be "
                            "larger than or equal to parameter 'min'");
  }
  if (minval != -1.) {
    if (minval == maxval and force.size() != 1) {
      throw std::domain_error(
          "TabulatedPotential parameter 'force' must contain 1 element");
    }
    if (force.empty()) {
      throw std::domain_error("TabulatedPotential parameter 'force' must "
                              "contain at least 1 element");
    }
    if (force.size() != energy.size()) {
      throw std::invalid_argument("TabulatedPotential parameter 'force' must "
                                  "have the same size as parameter 'energy'");
    }
    invstepsize = static_cast<double>(force.size() - 1) / (maxval - minval);
  } else {
    invstepsize = 0.;
  }
  force_tab = force;
  energy_tab = energy;
}
