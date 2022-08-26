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
 *  Implementation of \ref ljcos2.hpp
 */
#include "ljcos2.hpp"

#ifdef LJCOS2
#include "nonbonded_interaction_data.hpp"

#include <cmath>
#include <stdexcept>

LJcos2_Parameters::LJcos2_Parameters(double eps, double sig, double offset,
                                     double w)
    : eps{eps}, sig{sig}, offset{offset}, w{w} {
  if (eps < 0.) {
    throw std::domain_error("LJcos2 parameter 'epsilon' has to be >= 0");
  }
  if (sig < 0.) {
    throw std::domain_error("LJcos2 parameter 'sigma' has to be >= 0");
  }
  if (sig != 0.) {
    rchange = std::pow(2., 1. / 6.) * sig;
    cut = w + rchange;
  }
}

#endif /* ifdef LJCOS2 */
