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

#include <utils/math/sqr.hpp>

#include <numbers>
#include <stdexcept>

LJcos_Parameters::LJcos_Parameters(double epsilon, double sigma, double cutoff,
                                   double offset)
    : eps{epsilon}, sig{sigma}, cut{cutoff}, offset{offset} {
  if (epsilon < 0.) {
    throw std::domain_error("LJcos parameter 'epsilon' has to be >= 0");
  }
  if (sigma < 0.) {
    throw std::domain_error("LJcos parameter 'sigma' has to be >= 0");
  }
  if (cutoff < 0.) {
    throw std::domain_error("LJcos parameter 'cutoff' has to be >= 0");
  }

  constexpr auto cbrt2 = 1.25992104989487316476721060727822835057025;
  auto const facsq = cbrt2 * Utils::sqr(sig);

  rmin = sqrt(cbrt2) * sig;
  alfa = std::numbers::pi / (Utils::sqr(cut) - facsq);
  beta = std::numbers::pi * (1. - (1. / (Utils::sqr(cut) / facsq - 1.)));
}

#endif // LJCOS
