/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef UTILS_INTERPOLATION_DETAIL_LL_AND_DIST_HPP
#define UTILS_INTERPOLATION_DETAIL_LL_AND_DIST_HPP

#include "utils/Vector.hpp"

#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

namespace Utils {
namespace Interpolation {
namespace detail {

struct Block {
  /* Index of the lower left corner of the assignment cube */
  const std::array<int, 3> corner;
  /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
  const Vector3d distance;
};

template <size_t order, typename = void> struct ll_and_dist_;

template <size_t order>
struct ll_and_dist_<order, std::enable_if_t<(order % 2) == 1>> {
  Block operator()(const Vector3d &pos, const Vector3d &grid_spacing,
                   const Vector3d &offset) const {
    Vector3d dist;
    std::array<int, 3> ll;

    for (int dim = 0; dim < 3; dim++) {
      auto const fractional_index =
          (pos[dim] - offset[dim]) / grid_spacing[dim];

      const int nmp = std::floor(fractional_index + 0.5);
      dist[dim] = fractional_index - nmp;
      ll[dim] = nmp - (order - 1) / 2;
    }
    return {ll, dist};
  }
};

template <size_t order>
struct ll_and_dist_<order, std::enable_if_t<(order % 2) == 0>> {
  Block operator()(const Vector3d &pos, const Vector3d &grid_spacing,
                   const Vector3d &offset) const {
    Vector3d dist;
    std::array<int, 3> ll;

    for (int dim = 0; dim < 3; dim++) {
      auto const fractional_index =
          (pos[dim] - offset[dim]) / grid_spacing[dim];

      const int nmp = std::floor(fractional_index);
      dist[dim] = fractional_index - nmp - 0.5;
      ll[dim] = nmp - (order - 1) / 2;
    }
    return {ll, dist};
  }
};

/**
 * @brief Calculate the lower left index of a block
 *        stencil with order points side length.
 */
template <size_t order>
Block ll_and_dist(const Vector3d &pos, const Vector3d &grid_spacing,
                  const Vector3d &offset) {
  return ll_and_dist_<order>{}(pos, grid_spacing, offset);
}
} // namespace detail
} // namespace Interpolation
} // namespace Utils

#endif
