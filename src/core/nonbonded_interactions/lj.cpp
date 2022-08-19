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
 *  Implementation of \ref lj.hpp
 */
#include "lj.hpp"

#ifdef LENNARD_JONES
#include "nonbonded_interaction_data.hpp"

#include <utils/math/int_pow.hpp>

#include <algorithm>
#include <stdexcept>

LJ_Parameters::LJ_Parameters(double eps, double sig, double cut, double offset,
                             double min)
    : LJ_Parameters(eps, sig, cut, offset, min, 0.) {
  if (cut != 0.) {
    auto const sig_cut = sig / cut;
    shift = Utils::int_pow<6>(sig_cut) - Utils::int_pow<12>(sig_cut);
  }
}

LJ_Parameters::LJ_Parameters(double eps, double sig, double cut, double offset,
                             double min, double shift)
    : eps{eps}, sig{sig}, cut{cut}, shift{shift}, offset{offset}, min{std::max(
                                                                      min,
                                                                      0.)} {
  if (eps < 0.) {
    throw std::domain_error("LJ parameter 'epsilon' has to be >= 0");
  }
  if (sig < 0.) {
    throw std::domain_error("LJ parameter 'sigma' has to be >= 0");
  }
}

#endif /* ifdef LENNARD_JONES */
