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

#include "core/system/System.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameter.hpp>
#include <script_interface/system/Leaf.hpp>

#include <memory>
#include <optional>
#include <string>

namespace ScriptInterface::Coulomb {

class Container : public AutoParameters<Container, System::Leaf> {
  ObjectRef m_solver;
  ObjectRef m_extension;
  std::unique_ptr<VariantMap> m_params;

  void reset_solver() {
    auto &system = get_system();
    auto &coulomb = system.coulomb;
    m_solver.reset();
    m_extension.reset();
    coulomb.impl->extension = std::nullopt;
    coulomb.impl->solver = std::nullopt;
    system.on_coulomb_change();
  }

  void reset_extension() {
    auto &system = get_system();
    auto &coulomb = system.coulomb;
    m_extension.reset();
    coulomb.impl->extension = std::nullopt;
    system.on_coulomb_change();
  }

  void on_bind_system(::System::System &system) override {
    static_cast<void>(system);
    auto const &params = *m_params;
    for (auto const &key : get_parameter_insertion_order()) {
      if (params.count(key)) {
        do_set_parameter(key.c_str(), params.at(key));
      }
    }
    m_params.reset();
  }

public:
  Container() {
    add_parameters({
        {"solver",
         [this](Variant const &v) {
           if (m_extension) {
             if (context()->is_head_node()) {
               throw std::runtime_error(
                   "Cannot change solver when an extension is active");
             }
             throw Exception("");
           }
           if (is_none(v)) {
             reset_solver();
           } else {
             auto actor = get_value<ObjectRef>(v);
             auto leaf = std::dynamic_pointer_cast<System::Leaf>(actor);
             assert(leaf);
             leaf->bind_system(m_system.lock());
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
             auto leaf = std::dynamic_pointer_cast<System::Leaf>(actor);
             assert(leaf);
             leaf->bind_system(m_system.lock());
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
    m_params = std::make_unique<VariantMap>(params);
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
