/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "Integrator.hpp"

#include "core/integrate.hpp"

#include "script_interface/ScriptInterface.hpp"

#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace Integrators {

Variant Integrator::integrate(VariantMap const &params) {
  auto const steps = get_value<int>(params, "steps");
  auto const recalc_forces_flag = get_value<bool>(params, "recalc_forces");
  auto const reuse_forces_flag = get_value<bool>(params, "reuse_forces");

  context()->parallel_try_catch([&]() {
    if (steps < 0) {
      throw std::domain_error("Parameter 'steps' must be positive");
    }
    if (recalc_forces_flag and reuse_forces_flag) {
      throw std::invalid_argument(
          "cannot reuse old forces and recalculate forces");
    }
  });

  auto constexpr update_accumulators = true;
  auto const reuse_forces = reuse_forces_flag
                                ? INTEG_REUSE_FORCES_ALWAYS
                                : INTEG_REUSE_FORCES_CONDITIONALLY;
  return ::integrate_with_signal_handler(steps, reuse_forces,
                                         update_accumulators);
}

Variant Integrator::do_call_method(std::string const &name,
                                   VariantMap const &params) {
  if (name == "integrate") {
    return integrate(params);
  }
  return none;
}

} // namespace Integrators
} // namespace ScriptInterface
