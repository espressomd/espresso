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

#include <cassert>
#include <memory>

namespace System {

class System;

/**
 * @brief Abstract class that represents a component of the system.
 *
 * See @ref SystemClassDesign for more details.
 */
template <typename Class> class Leaf {
protected:
  std::weak_ptr<System> m_system;

  auto &get_system() {
    auto const ptr = m_system.lock();
    assert(ptr);
    return *ptr;
  }

  auto &get_system() const {
    auto const ptr = m_system.lock();
    assert(ptr);
    return *ptr;
  }

public:
  void bind_system(std::shared_ptr<System> const &system) {
    assert(system);
    assert(m_system.expired() or m_system.lock() == system);
    m_system = system;
  }

  void detach_system([[maybe_unused]] std::shared_ptr<System> const &system) {
    assert(system);
    assert(not m_system.expired());
    assert(system == m_system.lock());
    m_system.reset();
  }
};

} // namespace System
