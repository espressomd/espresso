/*
 * Copyright (C) 2021 The ESPResSo project
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

#include "LatticeWalberla.hpp"

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>

#include "walberla_utils.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <stdexcept>
#include <utility>

LatticeWalberla::LatticeWalberla(Utils::Vector3i const &grid_dimensions,
                                 Utils::Vector3i const &node_grid,
                                 unsigned int n_ghost_layers)
    : m_grid_dimensions{grid_dimensions}, m_n_ghost_layers{n_ghost_layers} {
  using walberla::real_t;
  using walberla::uint_c;

  for (int i : {0, 1, 2}) {
    if (m_grid_dimensions[i] % node_grid[i] != 0) {
      throw std::runtime_error(
          "LB grid dimensions and mpi node grid are not compatible.");
    }
  }

  auto constexpr lattice_constant = real_t{1};
  auto const cells_block = Utils::hadamard_division(grid_dimensions, node_grid);

  m_blocks = walberla::blockforest::createUniformBlockGrid(
      // number of blocks in each direction
      uint_c(node_grid[0]), uint_c(node_grid[1]), uint_c(node_grid[2]),
      // number of cells per block in each direction
      uint_c(cells_block[0]), uint_c(cells_block[1]), uint_c(cells_block[2]),
      lattice_constant,
      // number of cpus per direction
      uint_c(node_grid[0]), uint_c(node_grid[1]), uint_c(node_grid[2]),
      // periodicity
      true, true, true);
}

[[nodiscard]] std::pair<Utils::Vector3d, Utils::Vector3d>
LatticeWalberla::get_local_domain() const {
  using walberla::to_vector3d;
  // We only have one block per mpi rank
  assert(++(m_blocks->begin()) == m_blocks->end());

  auto const ab = m_blocks->begin()->getAABB();
  return {to_vector3d(ab.min()), to_vector3d(ab.max())};
}

[[nodiscard]] bool
LatticeWalberla::node_in_local_domain(const Utils::Vector3i &node) const {
  // Note: Lattice constant =1, cell centers offset by .5
  return ::walberla::get_block_and_cell(*this, node, false) != boost::none;
}
[[nodiscard]] bool
LatticeWalberla::node_in_local_halo(const Utils::Vector3i &node) const {
  return ::walberla::get_block_and_cell(*this, node, true) != boost::none;
}
[[nodiscard]] bool
LatticeWalberla::pos_in_local_domain(const Utils::Vector3d &pos) const {
  return ::walberla::get_block(*this, pos, false) != nullptr;
}
[[nodiscard]] bool
LatticeWalberla::pos_in_local_halo(const Utils::Vector3d &pos) const {
  return ::walberla::get_block(*this, pos, true) != nullptr;
}
