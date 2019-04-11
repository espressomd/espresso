/*
  Copyright (C) 2018 The ESPResSo project

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

#ifndef UTILS_UNIFORM_HPP
#define UTILS_UNIFORM_HPP

#include <cinttypes>
#include <limits>

namespace Utils {
/**
 * @brief Unifomly map unsigned integer to double.
 */
constexpr inline double uniform(uint64_t in) {
  auto constexpr const max = std::numeric_limits<uint64_t>::max();
  auto constexpr const fac = 1. / (max + 1.);

  return fac * in + 0.5 * fac;
}

} // namespace Utils

#endif // ESPRESSO_UNIFORM_HPP
