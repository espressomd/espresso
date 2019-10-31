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
#ifndef ESPRESSO_ABS_HPP
#define ESPRESSO_ABS_HPP

#include "utils/device_qualifier.hpp"

#ifndef __CUDACC__
#include <cmath>
#endif

namespace Utils {
/**
 * @brief Return the absolute value of x.
 */
inline DEVICE_QUALIFIER double abs(double x) { return fabs(x); }

/**
 * @brief Return the absolute value of x.
 */
inline DEVICE_QUALIFIER float abs(float x) { return fabsf(x); }
} // namespace Utils

#endif // ESPRESSO_DEVICE_MATH_HPP
