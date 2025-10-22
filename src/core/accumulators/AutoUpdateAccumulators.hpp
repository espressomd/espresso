/*
 * Copyright (C) 2016-2024 The ESPResSo project
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

#include "AccumulatorBase.hpp"

#include "system/Leaf.hpp"

namespace Accumulators {

class AutoUpdateAccumulators : public System::Leaf<AutoUpdateAccumulators> {
public:
  /**
   * @brief Update accumulators.
   *
   * Checks for all auto update accumulators if
   * they need to be updated and if so does.
   *
   */
  void operator()(boost::mpi::communicator const &comm, int steps);
  int next_update() const;
  bool contains(AccumulatorBase const *) const;
  void add(AccumulatorBase *);
  void remove(AccumulatorBase *);

private:
  struct AutoUpdateAccumulator {
    explicit AutoUpdateAccumulator(AccumulatorBase *acc)
        : frequency(acc->delta_N()), counter(1), acc(acc) {}
    int frequency;
    int counter;
    AccumulatorBase *acc;
  };

  std::vector<AutoUpdateAccumulator> m_accumulators;
};

} // namespace Accumulators
