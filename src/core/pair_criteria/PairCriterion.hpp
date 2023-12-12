/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef ESPRESSO_SRC_CORE_PAIR_CRITERIA_PAIR_CRITERION_HPP
#define ESPRESSO_SRC_CORE_PAIR_CRITERIA_PAIR_CRITERION_HPP

#include "Particle.hpp"
#include "particle_node.hpp"

namespace PairCriteria {
/**
 * @brief Criterion which returns a true/false value for a pair of particles.
 */
class PairCriterion {
public:
  /** @brief Make a decision based on two particles */
  virtual bool decide(Particle const &p1, Particle const &p2) const = 0;
  /**
   * @brief Make a decision based on particle ids.
   * This can only run on the head node outside of the integration loop.
   */
  bool decide(int id1, int id2) const {
    // Retrieve particle data
    auto const &p1 = get_particle_data(id1);
    auto const &p2 = get_particle_data(id2);
    const bool res = decide(p1, p2);
    return res;
  }

  virtual ~PairCriterion() = default;
};
} // namespace PairCriteria

#endif
