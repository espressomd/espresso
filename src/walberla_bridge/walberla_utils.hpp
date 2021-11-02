/*
 * Copyright (C) 2020 The ESPResSo project
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
#ifndef WALBERLA_UTILS_H
#define WALBERLA_UTILS_H

#include "blockforest/StructuredBlockForest.h"
#include "field/GhostLayerField.h"
#include "core/math/Vector3.h"
#include "core/math/Matrix3.h"

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>

#include <boost/optional.hpp>

#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>

namespace walberla {

// Vector conversion helpers
inline Utils::Vector3d to_vector3d(const Vector3<real_t> v) {
  return Utils::Vector3d{double_c(v[0]), double_c(v[1]), double_c(v[2])};
}
inline Vector3<real_t> to_vector3(const Utils::Vector3d v) {
  return Vector3<real_t>{real_c(v[0]), real_c(v[1]), real_c(v[2])};
}
inline Utils::Vector6d to_vector6d(const Matrix3<real_t> m) {
  return Utils::Vector6d{double_c(m[0]), double_c(m[3]), double_c(m[4]),
                         double_c(m[6]), double_c(m[7]), double_c(m[8])};
}
inline Utils::Vector3i to_vector3i(const std::array<int, 3> v) {
  return Utils::Vector3i{v[0], v[1], v[2]};
}

// Helpers to retrieve blocks and cells
struct BlockAndCell {
  IBlock *block;
  Cell cell;
};

template <typename PosVector, typename BlockForest>
IBlock *get_block_extended(const PosVector &pos, BlockForest const &blocks,
                           unsigned int n_ghost_layers) {
  for (auto b = blocks->begin(); b != blocks->end(); ++b) {
    if (b->getAABB()
            .getExtended(real_c(n_ghost_layers))
            .contains(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]))) {
      return &(*b);
    }
  }
  // Cell not in local blocks
  return {};
}

template <typename BlockForest>
boost::optional<BlockAndCell>
get_block_and_cell(const Utils::Vector3i &node, bool consider_ghost_layers,
                   BlockForest const &blocks, unsigned int n_ghost_layers) {
  // Get block and local cell
  Cell global_cell{uint_c(node[0]), uint_c(node[1]), uint_c(node[2])};
  auto block = blocks->getBlock(global_cell, 0);
  // Return if we don't have the cell
  if (consider_ghost_layers and !block) {
    // Try to find a block which has the cell as ghost layer
    block = get_block_extended(node, blocks, n_ghost_layers);
  }
  if (!block)
    return {boost::none};

  // Transform coords to block local
  Cell local_cell;
  blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);
  return {{block, local_cell}};
}

template <typename BlockForest>
IBlock *get_block(const Utils::Vector3d &pos, bool consider_ghost_layers,
                  BlockForest const &blocks, unsigned int n_ghost_layers) {
  // Get block
  auto block = blocks->getBlock(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]));
  if (consider_ghost_layers and !block) {
    block = get_block_extended(pos, blocks, n_ghost_layers);
  }
  return block;
}

template <typename Function>
void interpolate_bspline_at_pos(Utils::Vector3d pos, Function f) {
  Utils::Interpolation::bspline_3d<2>(
      pos, f, Utils::Vector3d{1.0, 1.0, 1.0}, // grid spacing
      Utils::Vector3d::broadcast(.5));        // offset
}

} // namespace walberla

#endif
