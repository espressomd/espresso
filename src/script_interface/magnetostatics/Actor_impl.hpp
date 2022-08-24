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

#ifdef DIPOLES

#include "Actor.hpp"

#include "core/magnetostatics/dipoles.hpp"
#include "core/magnetostatics/registration.hpp"

#include "script_interface/auto_parameters/AutoParameter.hpp"

namespace ScriptInterface {
namespace Dipoles {

template <class SIClass, class CoreClass>
Variant Actor<SIClass, CoreClass>::do_call_method(std::string const &name,
                                                  VariantMap const &params) {
  if (name == "activate") {
    context()->parallel_try_catch([&]() { ::Dipoles::add_actor(actor()); });
    return {};
  }
  if (name == "deactivate") {
    context()->parallel_try_catch([&]() { ::Dipoles::remove_actor(actor()); });
    return {};
  }
  return {};
}

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DIPOLES
