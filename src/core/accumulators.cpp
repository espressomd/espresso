/*
  Copyright (C) 2016-2018 The ESPResSo project

  This file is part of ESPResSo.


  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "accumulators.hpp"

#include <boost/range/algorithm/remove_if.hpp>
#include <boost/range/numeric.hpp>

#include <cassert>
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

void auto_update(int steps) {
  for (auto &acc : auto_update_accumulators) {
    assert(steps <= acc.frequency);
    acc.counter -= steps;
    if (acc.counter <= 0) {
      acc.acc->update();
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

void auto_update_add(AccumulatorBase *acc) {
  assert(acc);
  auto_update_accumulators.emplace_back(acc);
}
void auto_update_remove(AccumulatorBase *acc) {
  auto_update_accumulators.erase(
      boost::remove_if(
          auto_update_accumulators,
          [acc](AutoUpdateAccumulator const &au) { return au.acc == acc; }),
      auto_update_accumulators.end());
}

} // namespace Accumulators
