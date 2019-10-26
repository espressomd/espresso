/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef UTILS_MATH_GAUSSIAN_HPP
#define UTILS_MATH_GAUSSIAN_HPP

#include <utils/Vector.hpp>

#include <cmath>

namespace Utils {
inline double gaussian(Vector3d x, Vector3d x0, double sigma) {
  return std::exp(-((x - x0).norm2() / (2. * sigma * sigma)));
}

inline Utils::Vector3d del_gaussian(Vector3d x, Vector3d x0, double sigma) {
  return -(x - x0) * gaussian(x, x0, sigma) / (sigma * sigma);
}
} // namespace Utils

#endif // UTILS_MATH_GAUSSIAN_HPP
