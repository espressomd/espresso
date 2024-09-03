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

#include "CollisionPair.hpp"

#include "Particle.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"

#include <utility>
#include <vector>

namespace CollisionDetection {

class Off {
public:
  Off() = default;

  auto cutoff() const { return INACTIVE_CUTOFF; }

  void initialize(System::System &) {}

  void handle_collisions(System::System &, std::vector<CollisionPair> &) {}
  bool detect_collision(Particle const &, Particle const &,
                        double const) const {
    return false;
  }
};

} // namespace CollisionDetection
#endif // COLLISION_DETECTION
