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
#ifndef SCRIPT_INTERFACE_LEES_EDWARDS_OSCILLATORY_SHEAR_HPP
#define SCRIPT_INTERFACE_LEES_EDWARDS_OSCILLATORY_SHEAR_HPP

#include "core/lees_edwards/lees_edwards.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <boost/variant.hpp>

#include <memory>

namespace ScriptInterface {
namespace LeesEdwards {

class OscillatoryShear : public Protocol {
public:
  OscillatoryShear()
      : m_protocol{std::make_shared<::LeesEdwards::ActiveProtocol>(
            ::LeesEdwards::OscillatoryShear())} {
    add_parameters(
        {{"amplitude",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_amplitude},
         {"omega",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_omega},
         {"time_0",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_time_0}});
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
};

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
