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

#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <domain_decomposition/IBlock.h>

#include "LatticeWalberla.hpp"

#include <boost/optional.hpp>

#include <memory>

namespace walberla {
// Helpers to retrieve blocks and cells
struct BlockAndCell {
  IBlock *block;
  Cell cell;
};

template <typename T>
IBlock *get_block_extended(std::shared_ptr<StructuredBlockForest> const &blocks,
                           Utils::Vector<T, 3> const &pos,
                           unsigned int n_ghost_layers) {
  for (auto block = blocks->begin(); block != blocks->end(); ++block) {
    if (block->getAABB()
            .getExtended(real_c(n_ghost_layers))
            .contains(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]))) {
      return &(*block);
    }
  }
  // Cell not in local blocks
  return nullptr;
}

inline boost::optional<BlockAndCell>
get_block_and_cell(::LatticeWalberla const &lattice,
                   Utils::Vector3i const &node, bool consider_ghost_layers) {
  // Get block and local cell
  auto const blocks = lattice.get_blocks();
  Cell global_cell{uint_c(node[0]), uint_c(node[1]), uint_c(node[2])};
  auto block = blocks->getBlock(global_cell, 0);
  // Return if we don't have the cell
  if (consider_ghost_layers and !block) {
    // Try to find a block which has the cell as ghost layer
    block = get_block_extended(blocks, node, lattice.get_ghost_layers());
  }
  if (!block)
    return {boost::none};

  // Transform coords to block local
  Cell local_cell;
  blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);
  return {{block, local_cell}};
}

inline IBlock *get_block(::LatticeWalberla const &lattice,
                         Utils::Vector3d const &pos,
                         bool consider_ghost_layers) {
  // Get block
  auto const blocks = lattice.get_blocks();
  auto block = blocks->getBlock(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]));
  if (consider_ghost_layers and !block) {
    block = get_block_extended(blocks, pos, lattice.get_ghost_layers());
  }
  return block;
}

} // namespace walberla
