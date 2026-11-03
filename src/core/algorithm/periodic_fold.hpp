/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <cmath>
#include <concepts>
#include <limits>
#include <utility>

namespace Algorithm {
/**
 * @brief Fold value into primary interval.
 *
 * @param x Value to fold
 * @param i Image count before folding
 * @param l Length of primary interval
 * @return x folded into [0, l) and number of folds.
 */
inline auto periodic_fold(std::floating_point auto x, std::integral auto i,
                          std::floating_point auto l) noexcept {
  static_assert(std::is_same_v<decltype(x), decltype(l)>);
  using limits = std::numeric_limits<decltype(i)>;
  using value_type = decltype(x);

  while ((x < value_type{0}) && (i > limits::min())) {
    x += l;
    --i;
  }

  while ((x >= l) && (i < limits::max())) {
    x -= l;
    ++i;
  }

  return std::make_pair(x, i);
}

/**
 * @brief Fold value into primary interval.
 *
 * @param x Value to fold
 * @param l Length of primary interval
 * @return x folded into [0, l).
 */
inline auto periodic_fold(std::floating_point auto x,
                          std::floating_point auto l) noexcept {
  using value_type = decltype(x);
#ifndef __FAST_MATH__
  /* Can't fold if either x or l is nan or inf. */
  if (std::isnan(x) or std::isnan(l) or std::isinf(x) or (l == value_type{0})) {
    return std::nan("");
  }
  if (std::isinf(l)) {
    return x;
  }
#endif // __FAST_MATH__

  while (x < value_type{0}) {
    x += l;
  }

  while (x >= l) {
    x -= l;
  }

  return x;
}
} // namespace Algorithm
