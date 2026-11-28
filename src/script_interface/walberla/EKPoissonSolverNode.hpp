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
#include "LatticeIndices.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include <utils/Vector.hpp>

#include <cassert>
#include <memory>
#include <string>

namespace ScriptInterface::walberla {

class EKPoissonSolverNode
    : public AutoParameters<EKPoissonSolverNode, LatticeIndices> {
  std::shared_ptr<::walberla::PoissonSolver> m_ek_poisson_solver;
  Utils::Vector3i m_index;
  Utils::Vector3i m_grid_size;
  double m_conv_potential;

public:
  EKPoissonSolverNode() {
    add_parameters(
        {{"_index", AutoParameter::read_only, [this]() { return m_index; }}});
  }

  void do_construct(VariantMap const &params) override {
    auto const ek_sip =
        get_value<std::shared_ptr<EKPoissonSolver>>(params, "parent_sip");
    m_ek_poisson_solver = ek_sip->get_instance();
    assert(m_ek_poisson_solver);
    m_conv_potential = ek_sip->get_conversion_factor_potential();
    m_grid_size = m_ek_poisson_solver->get_lattice().get_grid_dimensions();
    m_index = get_mapped_index(get_value<Utils::Vector3i>(params, "index"),
                               m_grid_size);
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;
};
} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
