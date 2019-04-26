/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef UTILS_LINEAR_INTERPOLATION_HPP
#define UTILS_LINEAR_INTERPOLATION_HPP

#include <cassert>

namespace Utils {
/** Linear interpolation between two data points.
 *  @param[in] table   Tabulated values, equally-spaced along the x-axis
 *  @param[in] hi      %Distance on the x-axis between tabulated values
 *  @param[in] offset  Position on the x-axis of the first tabulated value
 *  @param[in] x       Position on the x-axis at which to interpolate the value
 *  @return Interpolated value on the y-axis at @p x.
 */
template <typename T, typename Container>
T linear_interpolation(Container const &table, T hi, T offset, T x) {
  auto const dind = (x - offset) * hi;
  auto const ind = static_cast<int>(dind);
  assert(ind <= dind);
  assert((ind >= 0) and (ind < table.size()));
  auto const dx = dind - ind;

  /* linear interpolation between data points */
  return table[ind] * (T{1} - dx) + table[ind + 1] * dx;
}
} // namespace Utils

#endif
