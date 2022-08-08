/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_GALILEI_COMFIXED_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_GALILEI_COMFIXED_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "core/forces.hpp"
#include "core/galilei/ComFixed.hpp"

#include <memory>
#include <vector>

namespace ScriptInterface {
namespace Galilei {

class ComFixed : public AutoParameters<ComFixed> {
  std::shared_ptr<::ComFixed> m_comfixed;

  void do_construct(VariantMap const &params) override {
    ::comfixed = m_comfixed = std::make_shared<::ComFixed>();
    for (auto const &p : params) {
      do_set_parameter(p.first, p.second);
    }
  }

public:
  ComFixed() {
    add_parameters({{"types",
                     [this](Variant const &v) {
                       auto const p_types = get_value<std::vector<int>>(v);
                       m_comfixed->set_fixed_types(p_types);
                     },
                     [this]() { return m_comfixed->get_fixed_types(); }}});
  }
};

} // namespace Galilei
} // namespace ScriptInterface
#endif
