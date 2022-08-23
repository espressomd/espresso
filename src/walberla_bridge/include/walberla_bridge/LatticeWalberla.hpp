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
#ifndef ESPRESSO_WALBERLA_LATTICE_HPP
#define ESPRESSO_WALBERLA_LATTICE_HPP

#include "utils/Vector.hpp"

#include <memory>
#include <utility>

namespace walberla::blockforest {
// forward declare
class StructuredBlockForest;
} // namespace walberla::blockforest

/** Class that runs and controls the BlockForest in walberla
 */
class LatticeWalberla {
public:
  using Lattice_T = walberla::blockforest::StructuredBlockForest;

private:
  Utils::Vector3i m_grid_dimensions;
  unsigned int m_n_ghost_layers;

  /** Block forest */
  std::shared_ptr<Lattice_T> m_blocks;

public:
  LatticeWalberla(Utils::Vector3i const &grid_dimensions,
                  Utils::Vector3i const &node_grid,
                  unsigned int n_ghost_layers);

  // Grid, domain, halo
  [[nodiscard]] auto get_ghost_layers() const { return m_n_ghost_layers; }
  [[nodiscard]] auto get_grid_dimensions() const { return m_grid_dimensions; }
  [[nodiscard]] auto get_blocks() const { return m_blocks; }
  [[nodiscard]] std::pair<Utils::Vector3d, Utils::Vector3d>
  get_local_domain() const;

  [[nodiscard]] bool node_in_local_domain(const Utils::Vector3i &node) const;
  [[nodiscard]] bool node_in_local_halo(const Utils::Vector3i &node) const;
  [[nodiscard]] bool pos_in_local_domain(const Utils::Vector3d &pos) const;
  [[nodiscard]] bool pos_in_local_halo(const Utils::Vector3d &pos) const;
};

#endif // ESPRESSO_WALBERLA_LATTICE_HPP
