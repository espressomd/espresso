/*
 * Copyright (C) 2019-2020 The ESPResSo project
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

#include "../src/lattice_boltzmann/LBWalberlaImpl.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <functional>
#include <memory>
#include <vector>

class LBTestParameters {
public:
  unsigned int seed;
  double kT;
  double viscosity;
  double density;
  Utils::Vector3d box_dimensions;
  Utils::Vector3i grid_dimensions;
  std::shared_ptr<LatticeWalberla> lattice;
};

using LbGeneratorVector = std::vector<
    std::function<std::shared_ptr<LBWalberlaBase>(LBTestParameters const &)>>;

LbGeneratorVector unthermalized_lbs() {
  LbGeneratorVector lbs;

  // Unthermalized D3Q19 MRT
  lbs.push_back([](LBTestParameters const &params) {
    auto ptr = std::make_shared<walberla::LBWalberlaImpl<>>(
        params.lattice, params.viscosity, params.density);
    ptr->set_collision_model(0.0, params.seed);
    ptr->ghost_communication();
    return ptr;
  });
  return lbs;
}

LbGeneratorVector thermalized_lbs() {
  LbGeneratorVector lbs;

  // Thermalized D3Q19 MRT
  lbs.push_back([](LBTestParameters const &params) {
    auto ptr = std::make_shared<walberla::LBWalberlaImpl<>>(
        params.lattice, params.viscosity, params.density);
    ptr->set_collision_model(params.kT, params.seed);
    ptr->ghost_communication();
    return ptr;
  });
  return lbs;
}

LbGeneratorVector all_lbs() {
  auto lbs = unthermalized_lbs();
  auto thermalized = thermalized_lbs();
  lbs.insert(lbs.end(), thermalized.begin(), thermalized.end());
  return lbs;
}

// Disable printing of type which does not support it
BOOST_TEST_DONT_PRINT_LOG_VALUE(LbGeneratorVector::value_type)

std::vector<Utils::Vector3i>
all_nodes_incl_ghosts(LatticeWalberla const &lattice, bool with_ghosts = true) {
  auto const &grid_dimensions = lattice.get_grid_dimensions();
  auto const n_ghost_layers =
      (with_ghosts) ? static_cast<int>(lattice.get_ghost_layers()) : 0;
  std::vector<Utils::Vector3i> res;
  for (int x = -n_ghost_layers; x < grid_dimensions[0] + n_ghost_layers; x++) {
    for (int y = -n_ghost_layers; y < grid_dimensions[1] + n_ghost_layers;
         y++) {
      for (int z = -n_ghost_layers; z < grid_dimensions[2] + n_ghost_layers;
           z++) {
        res.push_back(Utils::Vector3i{x, y, z});
      }
    }
  }
  return res;
}

std::vector<Utils::Vector3i>
local_nodes_incl_ghosts(LatticeWalberla const &lattice,
                        bool with_ghosts = true) {
  auto const [left, right] = lattice.get_local_domain();
  auto const n_ghost_layers =
      (with_ghosts) ? static_cast<int>(lattice.get_ghost_layers()) : 0;
  std::vector<Utils::Vector3i> res;
  for (int x = static_cast<int>(left[0]) - n_ghost_layers;
       x < static_cast<int>(right[0]) + n_ghost_layers; x++) {
    for (int y = static_cast<int>(left[1]) - n_ghost_layers;
         y < static_cast<int>(right[1]) + n_ghost_layers; y++) {
      for (int z = static_cast<int>(left[2]) - n_ghost_layers;
           z < static_cast<int>(right[2]) + n_ghost_layers; z++) {
        res.push_back(Utils::Vector3i{x, y, z});
      }
    }
  }
  return res;
}

std::vector<Utils::Vector3i> corner_nodes(Utils::Vector3i const &n) {
  std::vector<Utils::Vector3i> res;
  for (int i : {0, n[0] - 1})
    for (int j : {0, n[1] - 1})
      for (int k : {0, n[2] - 1})
        res.emplace_back(Utils::Vector3i{i, j, k});
  return res;
}

#endif
