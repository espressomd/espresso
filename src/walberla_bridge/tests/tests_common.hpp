/*
 * Copyright (C) 2019-2023 The ESPResSo project
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

#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <utils/Vector.hpp>

#include <initializer_list>
#include <memory>
#include <vector>

struct LatticeTestParameters {
  Utils::Vector3d box_dimensions;
  Utils::Vector3i grid_dimensions;
  std::shared_ptr<LatticeWalberla> lattice;
};

inline auto all_nodes_incl_ghosts(LatticeWalberla const &lattice,
                                  bool with_ghosts = true) {
  auto const &grid_dimensions = lattice.get_grid_dimensions();
  auto const gl =
      (with_ghosts) ? static_cast<int>(lattice.get_ghost_layers()) : 0;
  std::vector<Utils::Vector3i> res;
  for (auto x = -gl; x < grid_dimensions[0] + gl; ++x) {
    for (auto y = -gl; y < grid_dimensions[1] + gl; ++y) {
      for (auto z = -gl; z < grid_dimensions[2] + gl; ++z) {
        res.push_back(Utils::Vector3i{x, y, z});
      }
    }
  }
  return res;
}

inline auto local_nodes_incl_ghosts(LatticeWalberla const &lattice,
                                    bool with_ghosts = true) {
  auto const [left, right] = lattice.get_local_grid_range();
  auto const gl =
      (with_ghosts) ? static_cast<int>(lattice.get_ghost_layers()) : 0;
  std::vector<Utils::Vector3i> res;
  for (auto x = left[0] - gl; x < right[0] + gl; ++x) {
    for (auto y = left[1] - gl; y < right[1] + gl; ++y) {
      for (auto z = left[2] - gl; z < right[2] + gl; ++z) {
        res.push_back(Utils::Vector3i{x, y, z});
      }
    }
  }
  return res;
}

inline auto corner_nodes(Utils::Vector3i const &n) {
  std::vector<Utils::Vector3i> res;
  for (auto i : {0, n[0] - 1}) {
    for (auto j : {0, n[1] - 1}) {
      for (auto k : {0, n[2] - 1}) {
        res.emplace_back(Utils::Vector3i{i, j, k});
      }
    }
  }
  return res;
}

#endif // WALBERLA
