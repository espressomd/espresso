/*
 * Copyright (C) 2013-2022 The ESPResSo project
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

#include "System.hpp"

#include "script_interface/ScriptInterface.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace System {

/**
 * @brief Facade for the real @ref System class.
 * All calls are forwarded to the real class.
 * This construct is needed to satisfy the requirements
 * of the checkpointing mechanism, which instantiates
 * the @c SystemFacade after all other script interface
 * objects have been instantiated.
 */
class SystemFacade : public ObjectHandle {
  std::shared_ptr<System> m_instance;

public:
  SystemFacade() {
    // create a dummy system to be able to read the list of parameters
    m_instance = std::make_shared<System>();
  }
  void do_construct(VariantMap const &) override{};
  Variant get_parameter(const std::string &name) const override {
    return m_instance->get_parameter(name);
  }
  void do_set_parameter(const std::string &name, const Variant &v) override {
    m_instance->do_set_parameter(name, v);
  }
  Utils::Span<const boost::string_ref> valid_parameters() const override {
    return m_instance->valid_parameters();
  }
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_system_handle") {
      return m_instance;
    }
    if (name == "set_system_handle") {
      m_instance = std::dynamic_pointer_cast<System>(
          get_value<ObjectRef>(params, "obj"));
      return {};
    }
    return m_instance->do_call_method(name, params);
  }
};

} // namespace System
} // namespace ScriptInterface
