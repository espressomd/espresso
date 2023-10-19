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

#ifdef ELECTROSTATICS

#include "Actor.hpp"

#include "core/actor/registration.hpp"
#include "core/electrostatics/coulomb.hpp"
#include "core/system/System.hpp"

#include "script_interface/auto_parameters/AutoParameter.hpp"

namespace ScriptInterface {
namespace Coulomb {

template <class SIClass, class CoreClass> Actor<SIClass, CoreClass>::Actor() {
  add_parameters({
      {"prefactor", AutoParameter::read_only,
       [this]() { return m_actor->prefactor; }},
      {"check_neutrality",
       [this](Variant const &value) {
         auto const flag = get_value<bool>(value);
         auto &tolerance = m_actor->charge_neutrality_tolerance;
         if (flag) {
           if (tolerance == -1.) {
             tolerance = m_actor->charge_neutrality_tolerance_default;
           }
         } else {
           tolerance = -1.;
         }
       },
       [this]() {
         auto const tolerance = m_actor->charge_neutrality_tolerance;
         return Variant{tolerance != -1.};
       }},
      {"charge_neutrality_tolerance",
       [this](Variant const &value) {
         auto &tolerance = m_actor->charge_neutrality_tolerance;
         if (is_none(value)) {
           tolerance = -1.;
         } else {
           auto const new_tolerance = get_value<double>(value);
           if (new_tolerance < 0.) {
             if (context()->is_head_node()) {
               throw std::domain_error(
                   "Parameter 'charge_neutrality_tolerance' must be >= 0");
             }
             throw Exception("");
           }
           tolerance = new_tolerance;
         }
       },
       [this]() {
         auto const tolerance = m_actor->charge_neutrality_tolerance;
         if (tolerance == -1.) {
           return make_variant(none);
         }
         return Variant{tolerance};
       }},
  });
}

template <class SIClass, class CoreClass>
Variant Actor<SIClass, CoreClass>::do_call_method(std::string const &name,
                                                  VariantMap const &params) {
  if (name == "activate") {
    context()->parallel_try_catch([&]() {
      auto &system = System::get_system();
      add_actor(context()->get_comm(), system.coulomb.impl->solver, m_actor,
                [&system]() { system.on_coulomb_change(); });
    });
    return {};
  }
  return {};
}

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ELECTROSTATICS
