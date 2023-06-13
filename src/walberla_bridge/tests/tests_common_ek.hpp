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

#include "../src/electrokinetics/EKinWalberlaImpl.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <functional>
#include <memory>
#include <vector>

struct EKTestParameters : public LatticeTestParameters {
  unsigned int seed;
  double kT;
  double density;
  double diffusion;
  double valency;
  bool advection;
  bool friction_coupling;
  Utils::Vector3d ext_efield;
  Utils::Vector3d box_dimensions;
  Utils::Vector3i grid_dimensions;
  std::shared_ptr<LatticeWalberla> lattice;
};

using EkGeneratorVector = std::vector<
    std::function<std::shared_ptr<EKinWalberlaBase>(EKTestParameters const &)>>;

inline EkGeneratorVector unthermalized_eks() {
  using EKImplementation = walberla::EKinWalberlaImpl<>;
  EkGeneratorVector eks;

  eks.push_back([](EKTestParameters const &params) {
    auto ptr = std::make_shared<EKImplementation>(
        params.lattice, params.diffusion, 0., params.valency, params.ext_efield,
        params.density, params.advection, params.friction_coupling);
    ptr->ghost_communication();
    return ptr;
  });
  return eks;
}

inline EkGeneratorVector thermalized_eks() {
  using EKImplementation = walberla::EKinWalberlaImpl<>;
  EkGeneratorVector eks;

  eks.push_back([](EKTestParameters const &params) {
    auto ptr = std::make_shared<EKImplementation>(
        params.lattice, params.diffusion, params.kT, params.valency,
        params.ext_efield, params.density, params.advection,
        params.friction_coupling);
    ptr->ghost_communication();
    return ptr;
  });
  return eks;
}

inline EkGeneratorVector all_eks() {
  auto eks = unthermalized_eks();
  auto thermalized = thermalized_eks();
  eks.insert(eks.end(), thermalized.begin(), thermalized.end());
  return eks;
}

// Disable printing of type which does not support it
BOOST_TEST_DONT_PRINT_LOG_VALUE(EkGeneratorVector::value_type)

#endif // WALBERLA
