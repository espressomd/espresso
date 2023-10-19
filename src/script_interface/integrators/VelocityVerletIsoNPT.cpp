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

#ifdef NPT

#include "VelocityVerletIsoNPT.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/integrate.hpp"
#include "core/npt.hpp"
#include "core/system/System.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

VelocityVerletIsoNPT::VelocityVerletIsoNPT() {
  add_parameters({
      {"ext_pressure", AutoParameter::read_only,
       [this]() { return get_instance().p_ext; }},
      {"piston", AutoParameter::read_only,
       [this]() { return get_instance().piston; }},
      {"direction", AutoParameter::read_only,
       [this]() { return get_instance().get_direction(); }},
      {"cubic_box", AutoParameter::read_only,
       [this]() { return get_instance().cubic_box; }},
  });
}

void VelocityVerletIsoNPT::do_construct(VariantMap const &params) {
  auto const ext_pressure = get_value<double>(params, "ext_pressure");
  auto const piston = get_value<double>(params, "piston");
  auto const cubic_box = get_value_or<bool>(params, "cubic_box", false);
  auto const direction = get_value_or<Utils::Vector3b>(
      params, "direction", Utils::Vector3b::broadcast(true));

  context()->parallel_try_catch([&]() {
    m_instance = std::make_shared<::NptIsoParameters>(ext_pressure, piston,
                                                      direction, cubic_box);
  });
}

void VelocityVerletIsoNPT::activate() const {
  ::nptiso = get_instance();
  set_integ_switch(INTEG_METHOD_NPT_ISO);
  System::get_system().on_thermostat_param_change();
}

} // namespace Integrators
} // namespace ScriptInterface

#endif // NPT
