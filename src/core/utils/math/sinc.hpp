/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017,2018 The ESPResSo
  project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef UTILS_MATH_SINC_HPP
#define UTILS_MATH_SINC_HPP

#include "utils/constants.hpp"
#include "utils/device_qualifier.hpp"
#include "utils/math/abs.hpp"

#include <cmath>

namespace Utils {
/**
 * @brief Calculates the sinc-function as sin(PI*x)/(PI*x).
 *
 * (same convention as in Hockney/Eastwood). In order to avoid
 * divisions by 0, arguments, whose modulus is smaller than epsi, will
 * be evaluated by an 8th order Taylor expansion of the sinc
 * function. Note that the difference between sinc(x) and this
 * expansion is smaller than 0.235e-12, if x is smaller than 0.1. (The
 * next term in the expansion is the 10th order contribution
 * PI^10/39916800 * x^10 = 0.2346...*x^12).  This expansion should
 * also save time, since it reduces the number of function calls to
 * sin().
 */
template <typename T> DEVICE_QUALIFIER T sinc(T d) {
  const constexpr T epsi = 0.1;

  const auto PId = pi<T>() * d;

  if (::Utils::abs(d) > epsi)
    return sin(PId) / PId;

  /** Coefficients of the Taylor expansion of sinc */
  const constexpr T c2 = -0.1666666666667e-0;
  const constexpr T c4 = 0.8333333333333e-2;
  const constexpr T c6 = -0.1984126984127e-3;
  const constexpr T c8 = 0.2755731922399e-5;

  const auto PId2 = PId * PId;
  return T(1) + PId2 * (c2 + PId2 * (c4 + PId2 * (c6 + PId2 * c8)));
}
} // namespace Utils

#endif
