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

#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>

#include "utils/types_conversion.hpp"

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <domain_decomposition/IBlock.h>

#include <utils/Vector.hpp>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

LatticeWalberla::LatticeWalberla(Utils::Vector3i const &grid_dimensions,
                                 Utils::Vector3i const &node_grid,
                                 Utils::Vector3i const &block_grid,
                                 unsigned int n_ghost_layers)
    : m_grid_dimensions{grid_dimensions}, m_node_grid{node_grid},
      m_n_ghost_layers{n_ghost_layers} {
  using walberla::real_t;
  using walberla::uint_c;

  for (auto const i : {0u, 1u, 2u}) {
    if (m_grid_dimensions[i] % node_grid[i] != 0) {
      throw std::runtime_error(
          "Lattice grid dimensions and MPI node grid are not compatible.");
    }
    if (m_grid_dimensions[i] % block_grid[i] != 0) {
      throw std::runtime_error(
          "Lattice grid dimensions and block grid are not compatible.");
    }
  }

  auto constexpr lattice_constant = real_t{1};
  auto const cells_per_block =
      Utils::hadamard_division(grid_dimensions, block_grid);

  m_blocks = walberla::blockforest::createUniformBlockGrid(
      // number of blocks in each direction
      uint_c(block_grid[0]), uint_c(block_grid[1]), uint_c(block_grid[2]),
      // number of cells per block in each direction
      uint_c(cells_per_block[0]), uint_c(cells_per_block[1]),
      uint_c(cells_per_block[2]), lattice_constant,
      // number of cpus per direction
      uint_c(node_grid[0]), uint_c(node_grid[1]), uint_c(node_grid[2]),
      // periodicity
      true, true, true,
      // keep global block information
      false);
  for (IBlock &block : *m_blocks) {
    m_cached_blocks.push_back(&block);
  }
}

[[nodiscard]] std::pair<Utils::Vector3d, Utils::Vector3d>
LatticeWalberla::get_local_domain() const {
  using walberla::to_vector3d;
  // Get upper and lower corner of BlockForest assigned to a mpi rank.
  // Since we can allocate multiple blocks per mpi rank,
  // the corners of all Blocks are compared.
  auto aa = to_vector3d(m_blocks->begin()->getAABB().min());
  auto bb = to_vector3d(m_blocks->begin()->getAABB().max());
  for (auto const &block : *m_blocks) {
    auto cc = block.getAABB();
    for (auto const i : {0u, 1u, 2u}) {
      aa[i] = std::min(aa[i], cc.min()[i]);
      bb[i] = std::max(bb[i], cc.max()[i]);
    }
  }
  return {aa, bb};
}

[[nodiscard]] std::pair<Utils::Vector3i, Utils::Vector3i>
LatticeWalberla::get_local_grid_range(bool with_halo) const {
  auto const [lower_corner, upper_corner] = get_local_domain();
  auto const gl = with_halo ? static_cast<int>(m_n_ghost_layers) : 0;
  auto const padding = Utils::Vector3i::broadcast(gl);
  return {walberla::convert_cell_corner_to_coord(lower_corner) - padding,
          walberla::convert_cell_corner_to_coord(upper_corner) + padding};
}

[[nodiscard]] Utils::Vector3i
LatticeWalberla::get_block_corner(IBlock const &block, bool lower) const {
  if (lower) {
    return walberla::get_min_corner(block);
  }
  return walberla::get_max_corner(block);
}

[[nodiscard]] bool
LatticeWalberla::node_in_local_domain(Utils::Vector3i const &node) const {
  // Note: Lattice constant =1, cell centers offset by .5
  return ::walberla::get_block_and_cell(*this, node, false) != std::nullopt;
}
[[nodiscard]] bool
LatticeWalberla::node_in_local_halo(Utils::Vector3i const &node) const {
  return ::walberla::get_block_and_cell(*this, node, true) != std::nullopt;
}
[[nodiscard]] bool
LatticeWalberla::pos_in_local_domain(Utils::Vector3d const &pos) const {
  return ::walberla::get_block(*this, pos, false) != nullptr;
}
[[nodiscard]] bool
LatticeWalberla::pos_in_local_halo(Utils::Vector3d const &pos) const {
  return ::walberla::get_block(*this, pos, true) != nullptr;
}

[[nodiscard]] Utils::Vector3i
LatticeWalberla::calc_grid_dimensions(Utils::Vector3d const &box_size,
                                      double agrid) {
  auto const grid_dimensions =
      Utils::Vector3i{{static_cast<int>(std::round(box_size[0] / agrid)),
                       static_cast<int>(std::round(box_size[1] / agrid)),
                       static_cast<int>(std::round(box_size[2] / agrid))}};
  for (auto const i : {0u, 1u, 2u}) {
    if (std::abs(grid_dimensions[i] * agrid - box_size[i]) / box_size[i] >
        std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error(
          "Box length not commensurate with agrid in direction " +
          std::to_string(i) + " length " + std::to_string(box_size[i]) +
          " agrid " + std::to_string(agrid));
    }
  }
  return grid_dimensions;
}
