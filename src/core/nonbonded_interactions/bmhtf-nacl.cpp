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
 *  Implementation of \ref bmhtf-nacl.hpp
 */
#include "bmhtf-nacl.hpp"

#ifdef BMHTF_NACL
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/math/int_pow.hpp>

#include <stdexcept>

BMHTF_Parameters::BMHTF_Parameters(double a, double b, double c, double d,
                                   double sig, double cutoff)
    : A{a}, B{b}, C{c}, D{d}, sig{sig}, cut{cutoff} {
  if (a < 0.) {
    throw std::domain_error("BMHTF parameter 'a' has to be >= 0");
  }
  if (c < 0.) {
    throw std::domain_error("BMHTF parameter 'c' has to be >= 0");
  }
  if (d < 0.) {
    throw std::domain_error("BMHTF parameter 'd' has to be >= 0");
  }
  computed_shift = C / Utils::int_pow<6>(cut) + D / Utils::int_pow<8>(cut) -
                   A * std::exp(B * (sig - cut));
}

#endif // BMHTF_NACL
