/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "LatticeModel.hpp"
#include "VTKHandle.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <walberla_bridge/electrokinetics/PoissonSolver.hpp>

#include <memory>
#include <unordered_map>

namespace ScriptInterface::walberla {

class EKPoissonVTKHandle : public VTKHandleBase<::walberla::PoissonSolver> {
  static std::unordered_map<std::string, int> const obs_map;

  std::unordered_map<std::string, int> const &get_obs_map() const override {
    return obs_map;
  }
};

class EKPoissonSolver
    : public LatticeModel<::walberla::PoissonSolver, EKPoissonVTKHandle> {
public:
  virtual std::shared_ptr<::walberla::PoissonSolver>
  get_instance() const noexcept = 0;

  ::LatticeModel::units_map
  get_lattice_to_md_units_conversion() const override {
    return {
        {"potential", 1.},
    };
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
