/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#pragma once

#include "utils/Vector.hpp"

#include "detail/ll_and_dist.hpp"
#include "utils/math/bspline.hpp"
#include "utils/math/tensor_product.hpp"

#include <array>
#include <cstddef>

namespace Utils {
namespace Interpolation {
/**
 * @brief cardinal B-spline gradient interpolation with internal iteration.
 *
 * @tparam order The number of interpolation points in each direction.
 * @tparam Kernel Callable that can be called at every point with the
 *         index of the point and its weights as arguments.
 *
 * @param pos The point at which the interpolation should be run.
 * @param kernel Function to run for every point.
 * @param grid_spacing The distance between the grid points.
 * @param offset Shift of the grid relative to the origin.
 */
template <int order, typename Kernel>
void bspline_3d_gradient(Vector3d const &pos, Kernel const &kernel,
                         Vector3d const &grid_spacing, Vector3d const &offset) {
  using Utils::bspline;

  /* The coordinates and relative distance of the assignment cube. */
  auto const block = detail::ll_and_dist<order>(pos, grid_spacing, offset);

  /* Precalc weights that are used multiple times. */
  std::array<double, order> w_y;
  std::array<double, order> w_z;
  std::array<double, order> dw_y;
  std::array<double, order> dw_z;
  for (int i = 0; i < order; i++) {
    w_y[i] = bspline<order>(i, block.distance[1]);
    w_z[i] = bspline<order>(i, block.distance[2]);
    dw_y[i] = bspline_d<order>(i, block.distance[1]) / grid_spacing[1];
    dw_z[i] = bspline_d<order>(i, block.distance[2]) / grid_spacing[2];
  }

  std::array<int, 3> ind;
  for (int i = 0; i < order; i++) {
    ind[0] = block.corner[0] + i;
    auto const w_x = bspline<order>(i, block.distance[0]);
    auto const dw_x = bspline_d<order>(i, block.distance[0]) / grid_spacing[0];
    for (int j = 0; j < order; j++) {
      ind[1] = block.corner[1] + j;
      for (int k = 0; k < order; k++) {
        ind[2] = block.corner[2] + k;
        kernel(ind, Vector3d{dw_x * w_y[j] * w_z[k], w_x * dw_y[j] * w_z[k],
                             w_x * w_y[j] * dw_z[k]});
      }
    }
  }
}

/**
 * @brief cardinal B-spline weighted sum.
 */
template <int order, typename T, typename Kernel>
T bspline_3d_gradient_accumulate(Vector3d const &pos, Kernel const &kernel,
                                 Vector3d const &grid_spacing,
                                 Vector3d const &offset, T const &init) {
  T value = init;
  bspline_3d_gradient<order>(
      pos,
      [&value, &kernel](const std::array<int, 3> &ind, const Vector3d &w) {
        value += Utils::tensor_product(kernel(ind), w);
      },
      grid_spacing, offset);

  return value;
}
} // namespace Interpolation
} // namespace Utils
