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
 *  Implementation of \ref ljgen.hpp
 */
#include "ljgen.hpp"

#ifdef LENNARD_JONES_GENERIC
#include "nonbonded_interaction_data.hpp"

#include <stdexcept>

LJGen_Parameters::LJGen_Parameters(double epsilon, double sigma, double cutoff,
                                   double shift, double offset,
#ifdef LJGEN_SOFTCORE
                                   double lam, double delta,
#endif
                                   double e1, double e2, double b1, double b2)
    : eps{epsilon}, sig{sigma}, cut{cutoff}, shift{shift}, offset{offset},
#ifdef LJGEN_SOFTCORE
      lambda{lam}, softrad{delta},
#endif
      a1{e1}, a2{e2}, b1{b1}, b2{b2} {
  if (epsilon < 0.) {
    throw std::domain_error("Generic LJ parameter 'epsilon' has to be >= 0");
  }
  if (sigma < 0.) {
    throw std::domain_error("Generic LJ parameter 'sigma' has to be >= 0");
  }
#ifdef LJGEN_SOFTCORE
  if (delta < 0.) {
    throw std::domain_error("Generic LJ parameter 'delta' has to be >= 0");
  }
  if (lam < 0. or lam > 1.) {
    throw std::domain_error(
        "Generic LJ parameter 'lam' has to be in the range [0, 1]");
  }
#endif // LJGEN_SOFTCORE
}

#endif // LENNARD_JONES_GENERIC
