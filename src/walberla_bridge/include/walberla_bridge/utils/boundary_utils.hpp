/*
 * Copyright (C) 2020-2023 The ESPResSo project
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

#include "walberla_utils.hpp"

#include <utils/Vector.hpp>

#include <boost/multi_array/multi_array_ref.hpp>

#include <cassert>
#include <cstddef>
#include <iterator>
#include <vector>

namespace walberla {

inline std::vector<Utils::Vector3d>
fill_3D_vector_array(std::vector<double> const &vec_flat,
                     Utils::Vector3i const &grid_size) {
  auto const n_grid_points =
      static_cast<std::size_t>(Utils::product(grid_size));
  assert(vec_flat.size() == 3u * n_grid_points or vec_flat.size() == 3u);
  std::vector<Utils::Vector3d> output_vector;
  output_vector.reserve(3u * n_grid_points);

  auto const vec_begin = std::begin(vec_flat);
  auto const vec_end = std::end(vec_flat);
  if (vec_flat.size() == 3u) {
    auto const uniform_vector = Utils::Vector3d(vec_begin, vec_end);
    output_vector.assign(n_grid_points, uniform_vector);
  } else {
    output_vector.reserve(n_grid_points);
    for (auto it = vec_begin; it < vec_end; it += 3u) {
      output_vector.emplace_back(Utils::Vector3d(it, it + 3u));
    }
  }

  return output_vector;
}

inline std::vector<double>
fill_3D_scalar_array(std::vector<double> const &vec_flat,
                     Utils::Vector3i const &grid_size) {
  auto const n_grid_points =
      static_cast<std::size_t>(Utils::product(grid_size));
  assert(vec_flat.size() == n_grid_points or vec_flat.size() == 1u);
  std::vector<double> output_vector;
  output_vector.reserve(n_grid_points);

  auto const vec_begin = std::begin(vec_flat);
  auto const vec_end = std::end(vec_flat);
  if (vec_flat.size() == 1u) {
    auto const uniform_value = vec_flat[0];
    output_vector.assign(n_grid_points, uniform_value);
  } else {
    output_vector.assign(vec_begin, vec_end);
  }

  return output_vector;
}

template <class BoundaryModel, class DataType>
void set_boundary_from_grid(BoundaryModel &boundary,
                            LatticeWalberla const &lattice,
                            std::vector<int> const &raster_flat,
                            std::vector<DataType> const &data_flat) {
  // reshape grids
  auto const grid_size = lattice.get_grid_dimensions();
  assert(raster_flat.size() == Utils::product(grid_size));
  boost::const_multi_array_ref<DataType, 3> data_grid(data_flat.data(),
                                                      grid_size);
  boost::const_multi_array_ref<int, 3> raster(raster_flat.data(), grid_size);

  auto const &blocks = lattice.get_blocks();
  for (auto block = blocks->begin(); block != blocks->end(); ++block) {
    auto const [size_i, size_j, size_k] = boundary.block_dims(*block);
    auto const offset = lattice.get_local_grid_range().first;
    auto const off_i = offset[0];
    auto const off_j = offset[1];
    auto const off_k = offset[2];
    // Get field data which knows about the indices
    // In the loop, x,y,z are in block-local coordinates
    auto const n_ghost_layers = lattice.get_ghost_layers();
    auto const ghosts = static_cast<int>(n_ghost_layers);
    for (int i = off_i - ghosts; i < size_i + off_i + ghosts; ++i) {
      for (int j = off_j - ghosts; j < size_j + off_j + ghosts; ++j) {
        for (int k = off_k - ghosts; k < size_k + off_k + ghosts; ++k) {
          auto const node = Utils::Vector3i{{i, j, k}};
          auto const idx = (node + grid_size) % grid_size;
          if (raster(idx)) {
            auto const bc = get_block_and_cell(lattice, node, true);
            boundary.set_node_value_at_boundary(node, data_grid(idx), *bc);
          }
        }
      }
    }
  }
}

} // namespace walberla
