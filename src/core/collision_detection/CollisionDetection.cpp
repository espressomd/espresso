/*
 * Copyright (C) 2011-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_COLLISION_DETECTION

#include "CollisionDetection.hpp"

#include "cell_system/CellStructure.hpp"
#include "system/System.hpp"

#include <memory>
#include <utility>
#include <variant>

namespace CollisionDetection {

void CollisionDetection::handle_collisions() {
  auto &system = get_system();
  if (m_protocol) {
    std::visit(
        [&system, this](auto &protocol) {
          protocol.handle_collisions(system, local_collision_queue);
        },
        *m_protocol);
  }
  clear_queue();
}

void CollisionDetection::initialize() {
  auto &system = get_system();
  if (m_protocol) {
    std::visit([&system](auto &protocol) { protocol.initialize(system); },
               *m_protocol);
  }
  system.on_short_range_ia_change();
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  system.cell_structure->clear_local_properties();
#endif
}

void CollisionDetection::set_protocol(
    std::shared_ptr<ActiveProtocol> protocol) {
  m_protocol = std::move(protocol);
  initialize();
}

void CollisionDetection::unset_protocol() {
  m_protocol = nullptr;
  auto &system = get_system();
  system.on_short_range_ia_change();
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  system.cell_structure->clear_local_properties();
#endif
}

} // namespace CollisionDetection

#endif // ESPRESSO_COLLISION_DETECTION
