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

#include <array>
#include <cmath>

namespace Utils::Interpolation::detail {

struct Block {
  /** Index of the lower left corner of the assignment cube */
  std::array<int, 3> corner;
  /** Distance to the nearest mesh point in units of agrid in [-0.5, 0.5) */
  std::array<double, 3> distance;
};

/**
 * @brief Calculate the lower left index of a block
 *        stencil with @c order points side length.
 */
template <int order>
auto ll_and_dist(Vector3d const &pos, Vector3d const &grid_spacing,
                 Vector3d const &offset) {
  Block block{};
  for (unsigned int dim = 0u; dim < 3u; ++dim) {
    auto const fractional_index = (pos[dim] - offset[dim]) / grid_spacing[dim];
    int nmp;
    if constexpr (order % 2 == 0) {
      nmp = static_cast<int>(std::floor(fractional_index));
      block.distance[dim] = fractional_index - nmp - 0.5;
    } else {
      nmp = static_cast<int>(std::floor(fractional_index + 0.5));
      block.distance[dim] = fractional_index - nmp;
    }
    block.corner[dim] = nmp - (order - 1) / 2;
  }
  return block;
}

} // namespace Utils::Interpolation::detail
