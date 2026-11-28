/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "EKPoissonSolver.hpp"
#include "EKSpeciesSlice.hpp"

#include "LatticeSlice.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <walberla_bridge/electrokinetics/PoissonSolver.hpp>

#include <utils/Vector.hpp>

#include <cassert>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {

class EKPoissonSolverSlice : public LatticeSlice<EKFieldSerializer> {
  using LatticeModel = ::walberla::PoissonSolver;
  std::shared_ptr<LatticeModel> m_ek_poisson_solver;
  std::shared_ptr<EKPoissonSolver> m_ek_solver_sip;
  double m_conv_potential;
  std::unordered_map<std::string, std::vector<int>> m_shape_val;

public:
  void do_construct(VariantMap const &params) override {
    m_ek_solver_sip =
        get_value<std::shared_ptr<EKPoissonSolver>>(params, "parent_sip");
    m_ek_poisson_solver = m_ek_solver_sip->get_instance();
    assert(m_ek_poisson_solver);
    m_conv_potential = m_ek_solver_sip->get_conversion_factor_potential();
    m_shape = get_value<std::vector<int>>(params, "shape");
    m_slice_lower_corner =
        get_value<Utils::Vector3i>(params, "slice_lower_corner");
    m_slice_upper_corner =
        get_value<Utils::Vector3i>(params, "slice_upper_corner");
    m_shape_val["potential"] = std::vector<int>{1};
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  ::LatticeWalberla const &get_lattice() const override {
    return m_ek_poisson_solver->get_lattice();
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
