/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include <algorithm>
#include <iterator>

namespace Utils {
/**
 * @brief Check whether a range contains a value.
 *
 * Re-implementation of <tt>std::ranges::contains()</tt> from C++23.
 *
 * @param rng The range to search in.
 * @param value The value to search for.
 *
 * @return True iff range contains the value.
 */
template <class Range, class T> bool contains(Range &&rng, T const &value) {
  using std::end;

  return std::ranges::find(rng, value) != end(rng);
}
} // namespace Utils
