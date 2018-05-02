/*
  Copyright (C) 2016,2017 The ESPResSo project

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
#include "integrate.hpp"

namespace Accumulators {
std::vector<std::shared_ptr<Accumulators::AccumulatorBase>>
    auto_update_accumulators;

void auto_update() {
  for (auto &c : auto_update_accumulators) {
    c->update();
  }
}
}
