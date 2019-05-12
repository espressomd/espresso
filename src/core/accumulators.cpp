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

#include <boost/range/algorithm/remove.hpp>

#include <vector>

namespace Accumulators {
std::vector<AccumulatorBase *>
    auto_update_accumulators;

void auto_update() {
  for (auto &c : auto_update_accumulators) {
    c->auto_update();
  }
}

bool auto_update_enabled() { return !auto_update_accumulators.empty(); }

void auto_update_add(AccumulatorBase *acc) {
  auto_update_accumulators.push_back(acc);
}
void auto_update_remove(AccumulatorBase *acc) {
  auto_update_accumulators.erase(
      boost::remove(auto_update_accumulators, acc),
      auto_update_accumulators.end()
  );
}

} // namespace Accumulators
