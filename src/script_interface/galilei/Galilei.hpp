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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_GALILEI_GALILEI_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_GALILEI_GALILEI_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "core/galilei/Galilei.hpp"

#include <memory>

namespace ScriptInterface {
namespace Galilei {

class Galilei : public AutoParameters<Galilei> {
  std::shared_ptr<::Galilei> m_galilei;

  void do_construct(VariantMap const &params) override {
    m_galilei = std::make_shared<::Galilei>();
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "kill_particle_motion") {
      auto const rotation = get_value_or<bool>(params, "rotation", false);
      m_galilei->kill_particle_motion(rotation);
    }
    if (name == "kill_particle_forces") {
      auto const torque = get_value_or<bool>(params, "torque", false);
      m_galilei->kill_particle_forces(torque);
    }
    if (name == "galilei_transform") {
      m_galilei->galilei_transform();
    }
    if (name == "system_CMS") {
      return m_galilei->calc_system_cms_position().as_vector();
    }
    if (name == "system_CMS_velocity") {
      return m_galilei->calc_system_cms_velocity().as_vector();
    }
    return {};
  }
};

} // namespace Galilei
} // namespace ScriptInterface

#endif
