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

#pragma once

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <memory>
#include <utility>

namespace walberla::blockforest {
// forward declare
class StructuredBlockForest;
} // namespace walberla::blockforest

/** Class that runs and controls the BlockForest in waLBerla. */
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
  [[nodiscard]] auto get_local_grid_range() const {
    auto const conversion = [](Utils::Vector3d const &pos) -> Utils::Vector3i {
      auto const dim =
          Utils::Vector3i{{static_cast<int>(pos[0]), static_cast<int>(pos[1]),
                           static_cast<int>(pos[2])}};
#ifndef NDEBUG
      for (auto const i : {0u, 1u, 2u}) {
        assert(std::abs(static_cast<double>(dim[i]) - pos[i]) < 1e-10);
      }
#endif
      return dim;
    };
    auto const [lower_corner, upper_corner] = get_local_domain();
    return std::make_pair(conversion(lower_corner), conversion(upper_corner));
  }

  [[nodiscard]] bool node_in_local_domain(Utils::Vector3i const &node) const;
  [[nodiscard]] bool node_in_local_halo(Utils::Vector3i const &node) const;
  [[nodiscard]] bool pos_in_local_domain(Utils::Vector3d const &pos) const;
  [[nodiscard]] bool pos_in_local_halo(Utils::Vector3d const &pos) const;
  [[nodiscard]] static Utils::Vector3i
  calc_grid_dimensions(Utils::Vector3d const &box_size, double agrid);
};
