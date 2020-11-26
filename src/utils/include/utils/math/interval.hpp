/*
 * Copyright (C) 2020 The ESPResSo project
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
#ifndef ESPRESSO_INTERVAL_HPP
#define ESPRESSO_INTERVAL_HPP

#include <cmath>

namespace Utils {
/**
 * @brief Wrap around the value of @p val in the interval <tt>[low, high]</tt>.
 */
template <typename T> T interval(T val, T low, T high) {
  if (val == low or val == high)
    return val;
  auto const span = high - low;
  return std::fmod(span + std::fmod(val - low, span), span) + low;
}
} // namespace Utils

#endif // ESPRESSO_INTERVAL_HPP
