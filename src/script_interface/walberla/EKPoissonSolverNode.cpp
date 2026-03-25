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

#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "EKPoissonSolverNode.hpp"

#include "LatticeIndices.hpp"

#include <walberla_bridge/electrokinetics/PoissonSolver.hpp>

#include <utils/Vector.hpp>
#include <utils/mpi/reduce_optional.hpp>

#include <string>

namespace ScriptInterface::walberla {

Variant EKPoissonSolverNode::do_call_method(std::string const &name,
                                            VariantMap const &params) {
  if (name == "override_index") {
    // this hidden feature is used to iterate an EK slice without
    // rebuilding an EKPoissonSolverNode for each node in the slice
    auto const index = get_value<Utils::Vector3i>(params, "index");
    if (not is_index_valid(index, m_grid_size)) {
      return 1;
    }
    m_index = index;
    return 0;
  }
  if (name == "get_potential") {
    auto const result = m_ek_poisson_solver->get_node_potential(m_index);
    return Utils::Mpi::reduce_optional(context()->get_comm(), result);
  }

  return {};
}

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
