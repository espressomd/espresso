/*
 * Copyright (C) 2010-2020 The ESPResSo project
 * Copyright (C) 2002-2010
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
#ifndef P3M_DATA_STRUCT_HPP
#define P3M_DATA_STRUCT_HPP

#include "config.hpp"

#if defined(P3M) || defined(DP3M)

#include "p3m-common.hpp"

struct p3m_data_struct_base {
  P3MParameters params;

  /** Spatial differential operator in k-space. We use an i*k differentiation.
   */
  std::array<std::vector<int>, 3> d_op;
  /** Force optimised influence function (k-space) */
  std::vector<double> g_force;
  /** Energy optimised influence function (k-space) */
  std::vector<double> g_energy;

  /** number of permutations in k_space */
  int ks_pnum;

  /** Calculate the Fourier transformed differential operator.
   *  Remark: This is done on the level of n-vectors and not k-vectors,
   *  i.e. the prefactor @f$ 2i\pi/L @f$ is missing!
   */
  void calc_differential_operator() {
    d_op = detail::calc_meshift(
        {params.mesh[0], params.mesh[1], params.mesh[2]}, true);
  }
};

#endif
#endif
