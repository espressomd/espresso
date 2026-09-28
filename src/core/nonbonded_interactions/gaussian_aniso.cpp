/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
 *  Implementation of \ref gaussian_aniso.hpp
 */

#include "gaussian_aniso.hpp"

#ifdef ESPRESSO_GAUSSIAN_ANISO
#include "nonbonded_interaction_data.hpp"

#include <stdexcept>

GaussianAniso_Parameters::GaussianAniso_Parameters(double eps, double sig_x,
                                                   double sig_y, double sig_z,
                                                   double cutoff)
    : eps{eps}, sig_x{sig_x}, sig_y{sig_y}, sig_z{sig_z}, cut{cutoff} {
  if (sig_x <= 0.) {
    throw std::domain_error("GaussianAniso parameter 'sig_x' has to be >= 0");
  }
  if (sig_y <= 0.) {
    throw std::domain_error("GaussianAniso parameter 'sig_y' has to be >= 0");
  }
  if (sig_z <= 0.) {
    throw std::domain_error("GaussianAniso parameter 'sig_z' has to be >= 0");
  }
  if (cutoff <= 0.) {
    throw std::domain_error("GaussianAniso parameter 'cutoff' has to be >= 0");
  }
}

#endif // ESPRESSO_GAUSSIAN_ANISO
