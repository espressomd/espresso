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
 *  Implementation of \ref ljcos.hpp
 */
#include "ljcos.hpp"

#ifdef LJCOS
#include "nonbonded_interaction_data.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <stdexcept>

LJcos_Parameters::LJcos_Parameters(double eps, double sig, double cut,
                                   double offset)
    : eps{eps}, sig{sig}, cut{cut}, offset{offset} {
  if (eps < 0.) {
    throw std::domain_error("LJcos parameter 'epsilon' has to be >= 0");
  }
  if (sig < 0.) {
    throw std::domain_error("LJcos parameter 'sigma' has to be >= 0");
  }
  auto const facsq = Utils::cbrt_2() * Utils::sqr(sig);

  rmin = sqrt(Utils::cbrt_2()) * sig;
  alfa = Utils::pi() / (Utils::sqr(cut) - facsq);
  beta = Utils::pi() * (1. - (1. / (Utils::sqr(cut) / facsq - 1.)));
}

#endif // LJCOS
