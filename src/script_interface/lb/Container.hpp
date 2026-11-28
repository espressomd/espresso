/*
 * Copyright (C) 2023-2026 The ESPResSo project
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

#include <config/config.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameter.hpp>
#include <script_interface/system/Leaf.hpp>

#include "core/system/System.hpp"

#include <memory>
#include <optional>
#include <string>

namespace ScriptInterface::LB {

class Container : public AutoParameters<Container, System::Leaf> {
  ObjectRef m_solver;
  std::unique_ptr<VariantMap> m_params;

  void detach_solver() {
    if (m_solver) {
      m_solver->do_call_method("deactivate", {});
      std::dynamic_pointer_cast<System::Leaf>(m_solver)->detach_system();
    }
    m_solver.reset();
  }

  void bind_solver() {
    auto system = m_system.lock();
    std::dynamic_pointer_cast<System::Leaf>(m_solver)->bind_system(system);
    m_solver->do_call_method("activate", {});
  }

  void on_bind_system(::System::System &) override {
    auto const &params = *m_params;
    for (auto const &key : get_parameter_insertion_order()) {
      if (params.contains(key)) {
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
           if (is_none(v)) {
             detach_solver();
           } else {
             auto new_solver = get_value<ObjectRef>(v);
             auto old_solver = m_solver;
             detach_solver();
             try {
               m_solver = new_solver;
               context()->parallel_try_catch([this]() { bind_solver(); });
             } catch (...) {
               detach_solver();
               m_solver = old_solver;
               if (m_solver) {
                 bind_solver();
               }
               throw;
             }
           }
         },
         [this]() { return m_solver ? Variant{m_solver} : Variant{None{}}; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

protected:
  Variant do_call_method(std::string const &name, VariantMap const &) override {
    if (name == "clear") {
      detach_solver();
      return {};
    }
    return {};
  }
};

} // namespace ScriptInterface::LB
