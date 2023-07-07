/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "core/event.hpp"
#include "core/system/System.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameter.hpp>

#include <memory>
#include <optional>
#include <string>

namespace ScriptInterface::Coulomb {

class Container : public AutoParameters<Container> {
  ObjectRef m_solver;
  ObjectRef m_extension;

  void reset_solver() {
    auto &coulomb = System::get_system().coulomb;
    m_solver.reset();
    m_extension.reset();
    coulomb.extension = std::nullopt;
    coulomb.solver = std::nullopt;
    ::on_coulomb_change();
  }

  void reset_extension() {
    auto &coulomb = System::get_system().coulomb;
    m_extension.reset();
    coulomb.extension = std::nullopt;
    ::on_coulomb_change();
  }

public:
  Container() {
    add_parameters({
        {"solver",
         [this](Variant const &v) {
           if (is_none(v)) {
             reset_solver();
           } else {
             auto actor = get_value<ObjectRef>(v);
             actor->do_call_method("activate", {});
             m_solver = actor;
           }
         },
         [this]() { return m_solver ? Variant{m_solver} : Variant{None{}}; }},
        {"extension",
         [this](Variant const &v) {
           if (is_none(v)) {
             reset_extension();
           } else {
             auto actor = get_value<ObjectRef>(v);
             actor->do_call_method("activate", {});
             m_extension = actor;
           }
         },
         [this]() {
           return m_extension ? Variant{m_extension} : Variant{None{}};
         }},
    });
  }

  void do_construct(VariantMap const &params) override {
    if (params.count("solver")) {
      do_set_parameter("solver", params.at("solver"));
    }
    if (params.count("extension")) {
      do_set_parameter("extension", params.at("extension"));
    }
  }

protected:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "clear") {
      reset_solver();
      return {};
    }
    return {};
  }
};

} // namespace ScriptInterface::Coulomb

#endif // ELECTROSTATICS
