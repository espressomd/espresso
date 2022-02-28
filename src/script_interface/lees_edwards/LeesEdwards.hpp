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
#ifndef SCRIPT_INTERFACE_LEES_EDWARDS_LEES_EDWARDS_HPP
#define SCRIPT_INTERFACE_LEES_EDWARDS_LEES_EDWARDS_HPP

#include "Protocol.hpp"

#include "core/grid.hpp"
#include "core/lees_edwards.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>

namespace ScriptInterface {
namespace LeesEdwards {

class LeesEdwards : public AutoParameters<LeesEdwards> {
public:
  LeesEdwards() : m_protocol{nullptr} {
    add_parameters(
        {{"protocol",
          [this](Variant const &value) {
            if (is_none(value)) {
              m_protocol = nullptr;
              ::LeesEdwards::unset_protocol();
              return;
            }
            m_protocol = get_value<std::shared_ptr<Protocol>>(value);
            ::LeesEdwards::set_protocol(m_protocol->protocol());
          },
          [this]() {
            if (m_protocol)
              return make_variant(m_protocol);
            return make_variant(none);
          }},
         {"shear_velocity", box_geo.lees_edwards_bc().shear_velocity},
         {"pos_offset", box_geo.lees_edwards_bc().pos_offset},
         {"shear_direction", box_geo.lees_edwards_bc().shear_direction},
         {"shear_plane_normal", box_geo.lees_edwards_bc().shear_plane_normal}});
  }

private:
  std::shared_ptr<Protocol> m_protocol;
};

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
