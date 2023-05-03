/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <boost/multi_array.hpp>

#include <vector>

namespace Shapes {
std::vector<int> Shape::rasterize(Utils::Vector3i const &grid_size,
                                  double grid_spacing,
                                  double grid_offset) const {
  boost::multi_array<int, 3> raster(grid_size);
  for (int i = 0; i < grid_size[0]; ++i) {
    for (int j = 0; j < grid_size[1]; ++j) {
      for (int k = 0; k < grid_size[2]; ++k) {
        auto const pos = Utils::Vector3d{{(i + grid_offset) * grid_spacing,
                                          (j + grid_offset) * grid_spacing,
                                          (k + grid_offset) * grid_spacing}};
        raster[i][j][k] = is_inside(pos);
      }
    }
  }
  return {raster.data(), raster.data() + raster.num_elements()};
}
} // namespace Shapes
