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
/** \file
 *
 *  Implementation of \ref hat.hpp
 */
#include "hat.hpp"

#ifdef HAT
#include "nonbonded_interaction_data.hpp"

#include <stdexcept>

Hat_Parameters::Hat_Parameters(double F_max, double cutoff)
    : Fmax{F_max}, r{cutoff} {
  if (F_max < 0.) {
    throw std::domain_error("Hat parameter 'F_max' has to be >= 0");
  }
  if (cutoff < 0.) {
    throw std::domain_error("Hat parameter 'cutoff' has to be >= 0");
  }
}

#endif // HAT
