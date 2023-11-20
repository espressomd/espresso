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

/** @file
 *  Common parts of the MMM family of methods for the electrostatic
 *  interaction: MMM1D and ELC. This file contains the code for the
 *  polygamma expansions used for the near formulas of MMM1D.
 *
 *  The expansion of the polygamma functions is fairly easy and follows
 *  directly from @cite abramowitz65a. For details, see @cite arnold02a.
 */

#pragma once

#include "mmm-modpsi.hpp"

#include "specfunc.hpp"

#include <vector>

/** Modified polygamma for even order <tt>2*n, n >= 0</tt> */
inline double mod_psi_even(int n, double x) {
  return evaluateAsTaylorSeriesAt(modPsi[2 * n], x * x);
}

/** Modified polygamma for odd order <tt>2*n+1, n>= 0</tt> */
inline double mod_psi_odd(int n, double x) {
  return x * evaluateAsTaylorSeriesAt(modPsi[2 * n + 1], x * x);
}
