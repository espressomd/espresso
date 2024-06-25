/*
 * Copyright (C) 2016-2022 The ESPResSo project
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
#include "accumulators.hpp"

#include <boost/range/numeric.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

namespace Accumulators {
namespace {
struct AutoUpdateAccumulator {
  explicit AutoUpdateAccumulator(AccumulatorBase *acc)
      : frequency(acc->delta_N()), counter(1), acc(acc) {}
  int frequency;
  int counter;
  AccumulatorBase *acc;
};

std::vector<AutoUpdateAccumulator> auto_update_accumulators;
} // namespace

void auto_update(boost::mpi::communicator const &comm, int steps) {
  for (auto &acc : auto_update_accumulators) {
    assert(steps <= acc.frequency);
    acc.counter -= steps;
    if (acc.counter <= 0) {
      acc.acc->update(comm);
      acc.counter = acc.frequency;
    }

    assert(acc.counter > 0);
  }
}

int auto_update_next_update() {
  return boost::accumulate(auto_update_accumulators,
                           std::numeric_limits<int>::max(),
                           [](int a, AutoUpdateAccumulator const &acc) {
                             return std::min(a, acc.counter);
                           });
}

namespace detail {
struct MatchPredicate {
  AccumulatorBase const *m_acc;
  template <typename T> bool operator()(T const &a) const {
    return a.acc == m_acc;
  }
};
} // namespace detail

void auto_update_add(AccumulatorBase *acc) {
  assert(not auto_update_contains(acc));
  auto_update_accumulators.emplace_back(acc);
}

void auto_update_remove(AccumulatorBase *acc) {
  assert(auto_update_contains(acc));
  std::erase_if(auto_update_accumulators, detail::MatchPredicate{acc});
}

bool auto_update_contains(AccumulatorBase const *acc) noexcept {
  assert(acc);
  return std::ranges::any_of(auto_update_accumulators,
                             detail::MatchPredicate{acc});
}

} // namespace Accumulators
