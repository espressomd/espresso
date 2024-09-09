/*
 * Copyright (C) 2011-2024 The ESPResSo project
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

#ifdef COLLISION_DETECTION

#include "ActiveProtocol.hpp"
#include "CollisionPair.hpp"

#include "Particle.hpp"
#include "system/Leaf.hpp"

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace CollisionDetection {

class CollisionDetection : public System::Leaf<CollisionDetection> {
  std::shared_ptr<ActiveProtocol> m_protocol;

public:
  CollisionDetection() = default;
  /** @brief Get currently active collision protocol. */
  auto get_protocol() const { return m_protocol; }
  /** @brief Set a new collision protocol. */
  void set_protocol(std::shared_ptr<ActiveProtocol> protocol);
  /** @brief Delete the currently active collision protocol. */
  void unset_protocol();
  /** @brief Validate parameters and create particle types if needed. */
  void initialize();

  auto is_off() const {
    return m_protocol == nullptr or std::holds_alternative<Off>(*m_protocol);
  }

  auto cutoff() const {
    if (m_protocol == nullptr) {
      return INACTIVE_CUTOFF;
    }
    return std::visit([](auto const &protocol) { return protocol.cutoff(); },
                      *m_protocol);
  }

  /// Handle queued collisions
  void handle_collisions();

  void clear_queue() { local_collision_queue.clear(); }

  /** @brief Detect (and queue) a collision between the given particles. */
  void detect_collision(Particle const &p1, Particle const &p2,
                        double const dist_sq) {
    if (m_protocol) {
      bool collision_detected = std::visit(
          [&p1, &p2, dist_sq](auto const &protocol) {
            return protocol.detect_collision(p1, p2, dist_sq);
          },
          *m_protocol);

      if (collision_detected) {
        local_collision_queue.emplace_back(p1.id(), p2.id());
      }
    }
  }
  /// During force calculation, colliding particles are recorded in the queue.
  /// The queue is processed after force calculation, when it is safe to add
  /// particles.
  std::vector<CollisionPair> local_collision_queue;
};

} // namespace CollisionDetection

#endif // COLLISION_DETECTION
