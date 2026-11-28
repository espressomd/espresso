/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#include <utils/Vector.hpp>

#include <memory>
#include <utility>

// forward declarations
namespace walberla::blockforest {
class StructuredBlockForest;
} // namespace walberla::blockforest
namespace walberla::domain_decomposition {
class IBlock;
} // namespace walberla::domain_decomposition

/** Class that runs and controls the BlockForest in waLBerla. */
class LatticeWalberla {
public:
  using Lattice_T = walberla::blockforest::StructuredBlockForest;

private:
  Utils::Vector3i m_grid_dimensions;
  Utils::Vector3i m_node_grid;
  unsigned int m_n_ghost_layers;

  /** Block forest */
  std::shared_ptr<Lattice_T> m_blocks;
  using IBlock = walberla::domain_decomposition::IBlock;
  std::vector<IBlock *> m_cached_blocks;

public:
  LatticeWalberla(Utils::Vector3i const &grid_dimensions,
                  Utils::Vector3i const &node_grid,
                  Utils::Vector3i const &block_grid,
                  unsigned int n_ghost_layers);

  // Grid, domain, halo
  [[nodiscard]] auto get_ghost_layers() const { return m_n_ghost_layers; }
  [[nodiscard]] auto const &get_grid_dimensions() const {
    return m_grid_dimensions;
  }
  [[nodiscard]] auto const &get_node_grid() const { return m_node_grid; }
  [[nodiscard]] auto get_blocks() const { return m_blocks; }
  [[nodiscard]] auto const &get_cached_blocks() const {
    return m_cached_blocks;
  }
  [[nodiscard]] std::pair<Utils::Vector3d, Utils::Vector3d>
  get_local_domain() const;
  [[nodiscard]] std::pair<Utils::Vector3i, Utils::Vector3i>
  get_local_grid_range(bool with_halo = false) const;

  [[nodiscard]] Utils::Vector3i get_block_corner(IBlock const &block,
                                                 bool lower) const;
  [[nodiscard]] bool node_in_local_domain(Utils::Vector3i const &node) const;
  [[nodiscard]] bool node_in_local_halo(Utils::Vector3i const &node) const;
  [[nodiscard]] bool pos_in_local_domain(Utils::Vector3d const &pos) const;
  [[nodiscard]] bool pos_in_local_halo(Utils::Vector3d const &pos) const;
  [[nodiscard]] static Utils::Vector3i
  calc_grid_dimensions(Utils::Vector3d const &box_size, double agrid);
};
