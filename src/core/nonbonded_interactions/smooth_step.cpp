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
 *  Implementation of \ref smooth_step.hpp
 */
#include "smooth_step.hpp"

#ifdef SMOOTH_STEP
#include "nonbonded_interaction_data.hpp"

#include <stdexcept>

SmoothStep_Parameters::SmoothStep_Parameters(double eps, double sig,
                                             double cutoff, double d, int n,
                                             double k0)
    : eps{eps}, sig{sig}, cut{cutoff}, d{d}, n{n}, k0{k0} {
  if (eps < 0.) {
    throw std::domain_error("SmoothStep parameter 'eps' has to be >= 0");
  }
}

#endif // SMOOTH_STEP
