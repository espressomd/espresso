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

#include "config/config.hpp"

#ifdef STOKESIAN_DYNAMICS

#include "StokesianDynamics.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/PropagationMode.hpp"
#include "core/integrators/Propagation.hpp"
#include "core/stokesian_dynamics/sd_interface.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace Integrators {

StokesianDynamics::StokesianDynamics() {
  add_parameters({
      {"viscosity", AutoParameter::read_only,
       [this]() { return get_instance().viscosity; }},
      {"radii", AutoParameter::read_only,
       [this]() {
         return make_unordered_map_of_variants(get_instance().radii);
       }},
      {"lubrication", AutoParameter::read_only,
       [this]() {
         return static_cast<bool>(get_instance().flags &
                                  static_cast<int>(sd_flags::LUBRICATION));
       }},
      {"self_mobility", AutoParameter::read_only,
       [this]() {
         return static_cast<bool>(get_instance().flags &
                                  static_cast<int>(sd_flags::SELF_MOBILITY));
       }},
      {"pair_mobility", AutoParameter::read_only,
       [this]() {
         return static_cast<bool>(get_instance().flags &
                                  static_cast<int>(sd_flags::PAIR_MOBILITY));
       }},
      {"approximation_method", AutoParameter::read_only,
       [this]() {
         return std::string(
             (get_instance().flags & static_cast<int>(sd_flags::FTS)) ? "fts"
                                                                      : "ft");
       }},
  });
}

void StokesianDynamics::do_construct(VariantMap const &params) {
  context()->parallel_try_catch([&]() {
    int bitfield = 0;
    if (get_value_or<bool>(params, "self_mobility", true)) {
      bitfield |= static_cast<int>(sd_flags::SELF_MOBILITY);
    }
    if (get_value_or<bool>(params, "pair_mobility", true)) {
      bitfield |= static_cast<int>(sd_flags::PAIR_MOBILITY);
    }
    auto const approx =
        get_value_or<std::string>(params, "approximation_method", "fts");
    if (approx == "fts") {
      bitfield |= static_cast<int>(sd_flags::FTS);
    } else if (approx != "ft") {
      throw std::invalid_argument("Unknown approximation '" + approx + "'");
    }
    m_instance = std::make_shared<::StokesianDynamicsParameters>(
        get_value<double>(params, "viscosity"),
        get_value<std::unordered_map<int, double>>(params, "radii"), bitfield);
  });
}

void StokesianDynamics::activate() {
  context()->parallel_try_catch([&]() {
    register_integrator(get_instance());
    get_system().propagation->set_integ_switch(INTEG_METHOD_SD);
  });
}

} // namespace Integrators
} // namespace ScriptInterface

#endif // STOKESIAN_DYNAMICS
