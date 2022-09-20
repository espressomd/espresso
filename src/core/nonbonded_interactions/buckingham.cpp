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
 *  Implementation of \ref buckingham.hpp
 */
#include "buckingham.hpp"

#ifdef BUCKINGHAM
#include "nonbonded_interaction_data.hpp"

#include <stdexcept>

Buckingham_Parameters::Buckingham_Parameters(double a, double b, double c,
                                             double d, double cutoff,
                                             double discont, double shift)
    : A{a}, B{b}, C{c}, D{d}, cut{cutoff}, discont{discont}, shift{shift} {
  if (a < 0.) {
    throw std::domain_error("Buckingham parameter 'a' has to be >= 0");
  }
  if (b < 0.) {
    throw std::domain_error("Buckingham parameter 'b' has to be >= 0");
  }
  if (c < 0.) {
    throw std::domain_error("Buckingham parameter 'c' has to be >= 0");
  }
  if (d < 0.) {
    throw std::domain_error("Buckingham parameter 'd' has to be >= 0");
  }
  if (cutoff < 0.) {
    throw std::domain_error("Buckingham parameter 'cutoff' has to be >= 0");
  }

  /* Replace the Buckingham potential for interatomic distance less
     than or equal to discontinuity by a straight line (F1+F2*r) */
  auto const F = buck_force_r(A, B, C, D, discont);
  F1 = buck_energy_r(A, B, C, D, shift, discont) + discont * F;
  F2 = -F;
}

#endif // BUCKINGHAM
