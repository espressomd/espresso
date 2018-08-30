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
#ifndef UTILS_INTERPOLATION_DETAIL_POS_SHIFT_HPP
#define UTILS_INTERPOLATION_DETAIL_POS_SHIFT_HPP

#include <cmath>

namespace Utils {
namespace Interpolation {
namespace detail {
/**
 * @brief Calculates the shift of the aufpunkt relative to the mesh.
 *
 * For even meshes, distances are calculated relative to the middle
 * of two grid points, hence the distances are shifted by half a grid
 * constant. For odd size, distances are relative to the central point,
 * so there is no shift.
 */
template <unsigned order> constexpr double pos_shift() {
  return 0.5 * (order % 2);
}
} // namespace detail
} // namespace Interpolation
} // namespace Utils

#endif
