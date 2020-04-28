/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef UTILS_RASTER_HPP
#define UTILS_RASTER_HPP

#include "utils/Vector.hpp"
#include "utils/math/make_lin_space.hpp"

#include <boost/multi_array.hpp>

namespace Utils {

/**
 * @brief Raster a function over a regular 3d grid.
 *
 * This evaluates a function over a regular grid and
 * returns the function values at the grid points.
 *
 * @param offset Position of the lowest grid point
 * @param grid_spacing Grid constant
 * @param size Grid size
 * @param f Function to evaluate
 * @return Function values at the grid points.
 */
template <class T, class F>
auto raster(Vector<T, 3> const &offset, Vector<T, 3> const &grid_spacing,
            Vector3i size, F f) {
  using R = decltype(f((offset)));

  boost::multi_array<R, 3> res(size);

  auto const end = offset + Vector<T, 3>{grid_spacing[0] * size[0],
                                         grid_spacing[1] * size[1],
                                         grid_spacing[2] * size[2]};

  auto it = res.data();
  for (auto x : make_lin_space(offset[0], end[0], size[0], false))
    for (auto y : make_lin_space(offset[1], end[1], size[1], false))
      for (auto z : make_lin_space(offset[2], end[2], size[2], false)) {
        *it++ = f(Vector<T, 3>{x, y, z});
      }

  return res;
}
} // namespace Utils

#endif // UTILS_RASTER_HPP
