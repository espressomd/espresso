/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 * %Lattice Boltzmann D3Q19 model.
 */

#ifndef D3Q19_H
#define D3Q19_H

#include <array>
#include <cstddef>

namespace D3Q19 {

static constexpr std::size_t n_vel = 19;

/** Velocity sub-lattice of the D3Q19 model */
static constexpr const std::array<std::array<double, 3>, 19> c = {
    {{{0., 0., 0.}},
     {{1., 0., 0.}},
     {{-1., 0., 0.}},
     {{0., 1., 0.}},
     {{0., -1., 0.}},
     {{0., 0., 1.}},
     {{0., 0., -1.}},
     {{1., 1., 0.}},
     {{-1., -1., 0.}},
     {{1., -1., 0.}},
     {{-1., 1., 0.}},
     {{1., 0., 1.}},
     {{-1., 0., -1.}},
     {{1., 0., -1.}},
     {{-1., 0., 1.}},
     {{0., 1., 1.}},
     {{0., -1., -1.}},
     {{0., 1., -1.}},
     {{0., -1., 1.}}}};

/** Coefficients for pseudo-equilibrium distribution of the D3Q19 model */
static constexpr const std::array<std::array<double, 4>, 19> coefficients = {
    {{{1. / 3., 1., 3. / 2., -1. / 2.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 18., 1. / 6., 1. / 4., -1. / 12.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}},
     {{1. / 36., 1. / 12., 1. / 8., -1. / 24.}}}};

/** Coefficients in the functional for the equilibrium distribution */
static constexpr const std::array<double, 19> w = {
    {1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18.,
     1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
     1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.}};

/* the following values are the (weighted) lengths of the vectors */
static constexpr const std::array<double, 19> w_k = {
    {1.0, 1. / 3., 1. / 3., 1. / 3., 2. / 3., 4. / 9., 4. / 3., 1. / 9.,
     1. / 9., 1. / 9., 2. / 3., 2. / 3., 2. / 3., 2. / 9., 2. / 9., 2. / 9.,
     2.0, 4. / 9., 4. / 3.}};

template <typename T>
static constexpr const T c_sound_sq = static_cast<T>(1. / 3.);

} // namespace D3Q19

#undef GCC_EXTERN_STATEMENT

#endif /* D3Q19_H */
