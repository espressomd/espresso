/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "SteepestDescent.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/integrate.hpp"
#include "core/integrators/steepest_descent.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

Variant SteepestDescent::integrate(VariantMap const &params) {
  auto constexpr reuse_forces = INTEG_REUSE_FORCES_NEVER;
  auto constexpr update_accumulators = false;
  auto const steps = get_value<int>(params, "steps");
  context()->parallel_try_catch([&]() {
    if (steps < 0) {
      throw std::domain_error("Parameter 'steps' must be positive");
    }
  });
  return ::integrate_with_signal_handler(steps, reuse_forces,
                                         update_accumulators);
}

SteepestDescent::SteepestDescent() {
  add_parameters({
      {"f_max", AutoParameter::read_only,
       [this]() { return get_instance().f_max; }},
      {"gamma", AutoParameter::read_only,
       [this]() { return get_instance().gamma; }},
      {"max_displacement", AutoParameter::read_only,
       [this]() { return get_instance().max_displacement; }},
  });
}

void SteepestDescent::do_construct(VariantMap const &params) {
  auto const f_max = get_value<double>(params, "f_max");
  auto const gamma = get_value<double>(params, "gamma");
  auto const max_d = get_value<double>(params, "max_displacement");

  context()->parallel_try_catch([&]() {
    m_instance =
        std::make_shared<::SteepestDescentParameters>(f_max, gamma, max_d);
  });
}

void SteepestDescent::activate() const {
  register_integrator(get_instance());
  set_integ_switch(INTEG_METHOD_STEEPEST_DESCENT);
}

} // namespace Integrators
} // namespace ScriptInterface
