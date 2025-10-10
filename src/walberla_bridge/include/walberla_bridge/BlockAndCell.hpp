/*
 * Copyright (C) 2020-2025 The ESPResSo project
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

#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <domain_decomposition/IBlock.h>
#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include "../src/utils/types_conversion.hpp"
#include "LatticeWalberla.hpp"

#include <array>
#include <cmath>
#include <concepts>
#include <memory>
#include <optional>
#include <type_traits>

namespace detail {
template <typename T> struct is_real_vector : std::false_type {};

template <std::floating_point T>
struct is_real_vector<std::array<T, 3>> : std::true_type {};

template <std::floating_point T>
struct is_real_vector<walberla::Vector3<T>> : std::true_type {};

template <std::floating_point T>
struct is_real_vector<Utils::Vector<T, 3>> : std::true_type {};
} // namespace detail

template <typename T>
concept real_vector = detail::is_real_vector<T>::value;

namespace detail {
template <typename T> struct is_signed_integral_vector : std::false_type {};

template <std::signed_integral T>
struct is_signed_integral_vector<std::array<T, 3>> : std::true_type {};

template <std::signed_integral T>
struct is_signed_integral_vector<walberla::Vector3<T>> : std::true_type {};

template <std::signed_integral T>
struct is_signed_integral_vector<Utils::Vector<T, 3>> : std::true_type {};

template <> struct is_signed_integral_vector<walberla::Cell> : std::true_type {
  static_assert(std::integral<walberla::cell_idx_t> and
                std::is_signed_v<walberla::cell_idx_t>);
};
} // namespace detail

template <typename T>
concept signed_integral_vector = detail::is_signed_integral_vector<T>::value;

namespace walberla {

inline Cell to_cell(signed_integral_vector auto const &xyz) {
  return {xyz[0], xyz[1], xyz[2]};
}

struct BlockAndCell {
  IBlock *block;
  Cell cell;
};

IBlock *get_block_extended(LatticeWalberla const &lattice, auto const &pos,
                           unsigned int n_ghost_layers) {
  auto const &cached_blocks = lattice.get_cached_blocks();
  for (auto &block : cached_blocks) {
    if (block->getAABB()
            .getExtended(real_c(n_ghost_layers))
            .contains(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]))) {
      return &(*block);
    }
  }
  // Cell not in local blocks
  return nullptr;
}

inline std::optional<BlockAndCell>
get_block_and_cell(::LatticeWalberla const &lattice,
                   signed_integral_vector auto const &node,
                   bool consider_ghost_layers) {
  auto const &blocks = lattice.get_blocks();
  auto n_ghost_layers = 0u;
  if (consider_ghost_layers) {
    n_ghost_layers = lattice.get_ghost_layers();
  }

  auto block = get_block_extended(lattice, node, n_ghost_layers);
  if (!block)
    return std::nullopt;

  // Transform coords to block local
  Cell local_cell;

  Cell global_cell = to_cell(node);
  blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);
  return {{block, local_cell}};
}

inline IBlock *get_block(::LatticeWalberla const &lattice,
                         real_vector auto const &pos,
                         bool consider_ghost_layers) {
  // Get block
  auto const blocks = lattice.get_blocks();
  auto block = blocks->getBlock(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]));
  if (consider_ghost_layers and !block) {
    block = get_block_extended(lattice, pos, lattice.get_ghost_layers());
  }
  return block;
}

/**
 * @brief Get the block-local coordinates of a block corner.
 *
 * This method leverages the fact that the grid spacing is unity in LB units,
 * i.e. floating-point coordinates can be cast to integers indices.
 */
inline auto convert_cell_corner_to_coord(real_vector auto const &corner) {
  return Utils::Vector3i{{static_cast<int>(std::round(corner[0])),
                          static_cast<int>(std::round(corner[1])),
                          static_cast<int>(std::round(corner[2]))}};
}

/** @brief Get the block-local coordinates of the lower corner of a block. */
inline auto get_min_corner(IBlock const &block) {
  return convert_cell_corner_to_coord(block.getAABB().minCorner());
}

/** @brief Get the block-local coordinates of the upper corner of a block. */
inline auto get_max_corner(IBlock const &block) {
  return convert_cell_corner_to_coord(block.getAABB().maxCorner());
}

[[nodiscard]] inline std::optional<walberla::cell::CellInterval>
get_interval(::LatticeWalberla const &lattice,
             Utils::Vector3i const &lower_corner,
             Utils::Vector3i const &upper_corner) {
  auto const &cell_min = lower_corner;
  auto const cell_max = upper_corner - Utils::Vector3i::broadcast(1);
  auto const lower_bc = get_block_and_cell(lattice, cell_min, true);
  auto const upper_bc = get_block_and_cell(lattice, cell_max, true);
  if (not lower_bc or not upper_bc) {
    return std::nullopt;
  }

  auto const block_extent =
      get_min_corner(*upper_bc->block) - get_min_corner(*lower_bc->block);
  auto const global_lower_cell = lower_bc->cell;
  auto const global_upper_cell = upper_bc->cell + to_cell(block_extent);
  return {CellInterval(global_lower_cell, global_upper_cell)};
}

// Interval within local block
[[nodiscard]] inline std::optional<walberla::cell::CellInterval>
get_block_interval(::LatticeWalberla const &lattice,
                   Utils::Vector3i const &lower_corner,
                   Utils::Vector3i const &upper_corner,
                   Utils::Vector3i const &block_offset, IBlock const &block) {
  auto block_lower_corner = lattice.get_block_corner(block, true);
  if (not(upper_corner > block_lower_corner)) {
    return std::nullopt;
  }
  for (uint_t f = 0u; f < 3u; ++f) {
    block_lower_corner[f] = std::max(block_lower_corner[f], lower_corner[f]);
  }
  auto block_upper_corner = lattice.get_block_corner(block, false);
  if (not(block_upper_corner > lower_corner)) {
    return std::nullopt;
  }
  for (uint_t f = 0u; f < 3u; ++f) {
    block_upper_corner[f] = std::min(block_upper_corner[f], upper_corner[f]);
  }
  block_upper_corner -= Utils::Vector3i::broadcast(1);
  auto const block_lower_cell = to_cell(block_lower_corner - block_offset);
  auto const block_upper_cell = to_cell(block_upper_corner - block_offset);
  return {CellInterval(block_lower_cell, block_upper_cell)};
}

/**
 * @brief Synchronize data between a sliced block and a container.
 *
 * Synchronize data between two data buffers representing sliced matrices
 * with different memory layouts. The kernel takes as argument an index
 * for the flattened data buffer containing the serialized block slice,
 * an index for the flattened I/O buffer, and a block-local node position.
 *
 * @param bci           Cell interval of the local block within a 3D slice
 * @param ci            Cell interval of the entire lattice within a 3D slice
 * @param block_offset  Origin of the local block
 * @param lower_corner  Lower corner of the 3D slice
 * @param kernel        Function to execute on the two data buffers
 */
template <typename Kernel>
void copy_block_buffer(CellInterval const &bci, CellInterval const &ci,
                       Utils::Vector3i const &block_offset,
                       Utils::Vector3i const &lower_corner, Kernel &&kernel) {
  auto const local_grid = to_vector3i(ci.max() - ci.min() + Cell(1, 1, 1));
  auto const block_grid = to_vector3i(bci.max() - bci.min() + Cell(1, 1, 1));
  auto const lower_cell = bci.min();
  auto const upper_cell = bci.max();
  // In the loop, x,y,z are in block coordinates
  // The field data given in the argument knows about BlockForest
  // lattice indices from lower_corner to upper_corner. It is converted
  // to block coordinates
  for (auto x = lower_cell.x(), i = 0; x <= upper_cell.x(); ++x, ++i) {
    for (auto y = lower_cell.y(), j = 0; y <= upper_cell.y(); ++y, ++j) {
      for (auto z = lower_cell.z(), k = 0; z <= upper_cell.z(); ++z, ++k) {
        auto const node = block_offset + Utils::Vector3i{{x, y, z}};
        auto const local_index = Utils::get_linear_index(
            node - lower_corner, local_grid, Utils::MemoryOrder::ROW_MAJOR);
        auto const block_index = Utils::get_linear_index(
            i, j, k, block_grid, Utils::MemoryOrder::ROW_MAJOR);
        kernel(static_cast<unsigned>(block_index),
               static_cast<unsigned>(local_index), node);
      }
    }
  }
}

} // namespace walberla
