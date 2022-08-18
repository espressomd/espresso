/*
 * Copyright (C) 2018-2022 The ESPResSo project
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
/** \file
 *
 *  Implementation of \ref wca.hpp
 */
#include "wca.hpp"

#ifdef WCA
#include "nonbonded_interaction_data.hpp"

#include <cmath>
#include <stdexcept>

WCA_Parameters::WCA_Parameters(double eps, double sig) : eps{eps}, sig{sig} {
  if (eps < 0.) {
    throw std::domain_error("WCA parameter 'epsilon' has to be >= 0");
  }
  if (sig < 0.) {
    throw std::domain_error("WCA parameter 'sigma' has to be >= 0");
  }
  if (sig != 0.) {
    cut = sig * std::pow(2., 1. / 6.);
  }
}

#endif /* ifdef WCA */
