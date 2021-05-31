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
#ifndef UTILS_CONSTANTS_HPP
#define UTILS_CONSTANTS_HPP

#include "device_qualifier.hpp"

#include <boost/math/constants/constants.hpp>

namespace Utils {

/*************************************************************/
/** \name Mathematical, physical and chemical constants.     */
/*************************************************************/
/**@{*/

/**
 * @brief Ratio of diameter and circumference of a circle.
 */
template <class T = double> DEVICE_QUALIFIER constexpr T pi() {
  return T(3.14159265358979323846264338328L);
}

/**
 * @brief One over square root of pi.
 */
template <class T = double> DEVICE_QUALIFIER constexpr T sqrt_pi_i() {
  return T(0.56418958354775627928034964498L);
}

/**
 * @brief Euler-Mascheroni constant.
 */
template <class T = double> DEVICE_QUALIFIER constexpr T gamma() {
  return T(0.57721566490153286060651209008L);
}

/**
 * @brief Natural logarithm of 2.
 */
template <class T = double> DEVICE_QUALIFIER constexpr T ln_2() {
  return T(0.6931471805599453094172321214581766L);
}

/**
 * @brief Square root of 2.
 */
template <class T = double> DEVICE_QUALIFIER constexpr T sqrt_2() {
  return T(1.4142135623730950488016887242096981L);
}

/**
 * @brief Cube root of 2.
 */
template <class T = double> DEVICE_QUALIFIER constexpr T cbrt_2() {
  return T(1.25992104989487316476721060727822835057025L);
}

/**@}*/

/// error code if no error occurred
#define ES_OK 0
/// error code if an error occurred
#define ES_ERROR 1

} // namespace Utils

#endif
