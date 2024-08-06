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

#pragma once

#include "config/config.hpp"

#ifdef DIPOLES

#include "Actor.hpp"

#include "core/actor/registration.hpp"
#include "core/magnetostatics/dipoles.hpp"

#include "script_interface/auto_parameters/AutoParameter.hpp"

#include <cassert>
#include <string>

namespace ScriptInterface {
namespace Dipoles {

template <class SIClass, class CoreClass>
Variant Actor<SIClass, CoreClass>::do_call_method(std::string const &name,
                                                  VariantMap const &params) {
  assert(params.empty());
  if (name == "activate") {
    context()->parallel_try_catch([this]() {
      auto &system = get_system();
      add_actor(context()->get_comm(), m_system.lock(),
                system.dipoles.impl->solver, m_actor,
                [&system]() { system.on_dipoles_change(); });
    });
    return {};
  }
  return {};
}

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DIPOLES
