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

#include "script_interface/ObjectHandle.hpp"

#include "core/system/System.hpp"

#include <cassert>
#include <memory>

namespace ScriptInterface {
namespace System {

class Leaf : public ObjectHandle {
  virtual void on_bind_system(::System::System &) {}
  virtual void on_detach_system(::System::System &) {}

protected:
  std::weak_ptr<::System::System> m_system;

  auto const &get_system() const {
    auto const ptr = m_system.lock();
    assert(ptr);
    return *ptr;
  }

  auto &get_system() {
    auto const ptr = m_system.lock();
    assert(ptr);
    return *ptr;
  }

public:
  void bind_system(std::shared_ptr<::System::System> const &system) {
    assert(m_system.expired());
    m_system = system;
    on_bind_system(*system);
  }

  void detach_system() {
    auto const ptr = m_system.lock();
    assert(ptr != nullptr);
    on_detach_system(*ptr);
    m_system.reset();
  }
};

} // namespace System
} // namespace ScriptInterface
