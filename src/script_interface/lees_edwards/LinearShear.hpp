/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include "core/lees_edwards/lees_edwards.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>
#include <variant>

namespace ScriptInterface {
namespace LeesEdwards {

class LinearShear : public Protocol {
  using CoreClass = ::LeesEdwards::LinearShear;

public:
  LinearShear()
      : m_protocol{
            std::make_shared<::LeesEdwards::ActiveProtocol>(CoreClass())} {
    add_parameters(
        {{"initial_pos_offset",
          std::get<CoreClass>(*m_protocol).m_initial_pos_offset},
         {"shear_velocity", std::get<CoreClass>(*m_protocol).m_shear_velocity},
         {"time_0", std::get<CoreClass>(*m_protocol).m_time_0}});
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
};

} // namespace LeesEdwards
} // namespace ScriptInterface
