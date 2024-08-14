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

#include "AutoUpdateAccumulators.hpp"

#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <vector>

namespace Accumulators {

void AutoUpdateAccumulators::operator()(boost::mpi::communicator const &comm,
                                        int steps) {
  for (auto &acc : m_accumulators) {
    assert(steps <= acc.frequency);
    acc.counter -= steps;
    if (acc.counter <= 0) {
      acc.acc->update(comm);
      acc.counter = acc.frequency;
    }

    assert(acc.counter > 0);
  }
}

int AutoUpdateAccumulators::next_update() const {
  return std::accumulate(m_accumulators.begin(), m_accumulators.end(),
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

void AutoUpdateAccumulators::add(AccumulatorBase *acc) {
  assert(not contains(acc));
  auto const *this_system = &get_system();
  if (acc->has_same_system_handle(nullptr)) {
    acc->override_system_handle(this_system);
  } else if (not acc->has_same_system_handle(this_system)) {
    throw std::runtime_error("This accumulator is bound to another system");
  }
  m_accumulators.emplace_back(acc);
}

void AutoUpdateAccumulators::remove(AccumulatorBase *acc) {
  assert(contains(acc));
  std::erase_if(m_accumulators, detail::MatchPredicate{acc});
}

bool AutoUpdateAccumulators::contains(AccumulatorBase const *acc) const {
  assert(acc);
  return std::ranges::any_of(m_accumulators, detail::MatchPredicate{acc});
}

} // namespace Accumulators
