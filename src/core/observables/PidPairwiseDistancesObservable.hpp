/*
 * Copyright (C) 2025-2026 The ESPResSo project
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
#ifndef OBSERVABLES_PIDPAIRWISEDISTANCESOBSERVABLE_HPP
#define OBSERVABLES_PIDPAIRWISEDISTANCESOBSERVABLE_HPP

#include "PidObservable.hpp"

#include <stdexcept>
#include <vector>

namespace Observables {

/** @brief Calculate pairwise distances between two sets of particles. */
class PidPairwiseDistancesObservable : public PidObservable {
  std::vector<int> const m_target_ids;

public:
  PidPairwiseDistancesObservable(std::vector<int> const &ids,
                                 std::vector<int> const &target_ids)
      : PidObservable(ids), m_target_ids(target_ids) {}
  std::vector<int> const &target_ids() const { return m_target_ids; }
};

} // Namespace Observables
#endif
