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
#ifndef ESPRESSO_ACCUMULATORS_HPP
#define ESPRESSO_ACCUMULATORS_HPP

#include "accumulators/AccumulatorBase.hpp"

namespace Accumulators {
/**
 * @brief Update accumulators.
 *
 * Checks for all auto update accumulators if
 * they need to be updated and if so does.
 *
 */
void auto_update(int steps);
int auto_update_next_update();
void auto_update_add(AccumulatorBase *);
void auto_update_remove(AccumulatorBase *);

} // namespace Accumulators

#endif // ESPRESSO_ACCUMULATORS_HPP
