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
#include "bonded_interactions/bonded_tab.hpp"

#include "errorhandling.hpp"

#include <utils/constants.hpp>

#include <cassert>
#include <vector>

TabulatedBond::TabulatedBond(double min, double max,
                             std::vector<double> const &energy,
                             std::vector<double> const &force) {
  assert(max >= min);
  assert((max == min) || force.size() > 1);
  assert(force.size() == energy.size());

  auto tab_pot = this->pot = std::make_shared<TabulatedPotential>();

  /* set table limits */
  tab_pot->minval = min;
  tab_pot->maxval = max;

  tab_pot->invstepsize = static_cast<double>(force.size() - 1) / (max - min);

  tab_pot->force_tab = force;
  tab_pot->energy_tab = energy;
}

TabulatedDistanceBond::TabulatedDistanceBond(double min, double max,
                                             std::vector<double> const &energy,
                                             std::vector<double> const &force)
    : TabulatedBond(min, max, energy, force) {
  /* set table limits */
  this->pot->minval = min;
  this->pot->maxval = max;
}

TabulatedAngleBond::TabulatedAngleBond(double min, double max,
                                       std::vector<double> const &energy,
                                       std::vector<double> const &force)
    : TabulatedBond(min, max, energy, force) {
  /* set table limits */
  this->pot->minval = 0.0;
  this->pot->maxval = Utils::pi() + ROUND_ERROR_PREC;
}

TabulatedDihedralBond::TabulatedDihedralBond(double min, double max,
                                             std::vector<double> const &energy,
                                             std::vector<double> const &force)
    : TabulatedBond(min, max, energy, force) {
  /* set table limits */
  this->pot->minval = 0.0;
  this->pot->maxval = 2.0 * Utils::pi() + ROUND_ERROR_PREC;
}
