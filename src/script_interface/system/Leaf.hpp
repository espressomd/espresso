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

/**
 * @brief Script interface wrapper for a component of the system class.
 *
 * This class manages a core leaf object. A leaf can exist without a system,
 * in which case the managed object cannot be fully initialized
 * and has nothing to act upon.
 * Binding a leaf to a system triggers an initialization of the managed object.
 * Detaching a leaf may trigger a memory deallocation of the managed object
 * resources, without deallocating the managed object itself.
 * This behavior can be leveraged to implement move semantics.
 * See @ref SystemClassDesign for more details.
 */
class Leaf : public ObjectHandle {
  /** @brief Callback triggered upon binding a leaf to a system. */
  virtual void on_bind_system(::System::System &) {}
  /** @brief Callback triggered upon detaching a leaf from a system. */
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
    assert(m_system.expired() or m_system.lock() == system);
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
