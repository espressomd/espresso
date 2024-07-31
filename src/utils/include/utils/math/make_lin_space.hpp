/*
 * Copyright (C) 2019-2022 The ESPResSo project
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
#include <cstddef>
#include <ranges>

namespace Utils {
/**
 * @brief Equally spaced values in interval
 *
 * Returns a range of equally spaced values in
 * the range @p start to @p stop, like @c numpy.linspace().
 *
 * @tparam T floating point type
 * @param start Start value of the interval
 * @param stop End value of the interval
 * @param number Number of partition points
 * @param endpoint If true, the last point is @p stop.
 * @return Range of equally-spaced values
 */
template <class T>
auto make_lin_space(T start, T stop, std::size_t number, bool endpoint = true) {
  auto const dx = (stop - start) / T(number - endpoint);

  return std::ranges::views::transform(
      std::views::iota(std::size_t{0u}, number),
      [dx, start](std::size_t i) { return start + T(i) * dx; });
}
} // namespace Utils
