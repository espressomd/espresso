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

} // namespace walberla
