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

#include "tests_common.hpp"

#include "../src/lattice_boltzmann/LBWalberlaImpl.hpp"

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <functional>
#include <memory>
#include <vector>

struct LBTestParameters : public LatticeTestParameters {
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

inline LbGeneratorVector unthermalized_lbs() {
  using LBImplementation = walberla::LBWalberlaImpl<double, lbmpy::Arch::CPU>;
  LbGeneratorVector lbs;

  // Unthermalized D3Q19 MRT
  lbs.push_back([](LBTestParameters const &params) {
    auto ptr = std::make_shared<LBImplementation>(
        params.lattice, params.viscosity, params.density);
    ptr->set_collision_model(0.0, params.seed);
    ptr->ghost_communication();
    return ptr;
  });
  return lbs;
}

inline LbGeneratorVector thermalized_lbs() {
  using LBImplementation = walberla::LBWalberlaImpl<double, lbmpy::Arch::CPU>;
  LbGeneratorVector lbs;

  // Thermalized D3Q19 MRT
  lbs.push_back([](LBTestParameters const &params) {
    auto ptr = std::make_shared<LBImplementation>(
        params.lattice, params.viscosity, params.density);
    ptr->set_collision_model(params.kT, params.seed);
    ptr->ghost_communication();
    return ptr;
  });
  return lbs;
}

inline LbGeneratorVector all_lbs() {
  auto lbs = unthermalized_lbs();
  auto thermalized = thermalized_lbs();
  lbs.insert(lbs.end(), thermalized.begin(), thermalized.end());
  return lbs;
}

// Disable printing of type which does not support it
BOOST_TEST_DONT_PRINT_LOG_VALUE(LbGeneratorVector::value_type)

#endif // WALBERLA
