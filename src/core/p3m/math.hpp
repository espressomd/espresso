/*
 * Copyright (C) 2010-2024 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <utils/device_qualifier.hpp>
#include <utils/math/sqr.hpp>

#ifndef __CUDACC__
#include <cmath>
#endif
#include <numbers>
#include <stdexcept>
#include <string>

namespace math {

/** @brief Return the absolute value of x. */
inline DEVICE_QUALIFIER auto abs(double x) { return fabs(x); }

/** @brief Return the absolute value of x. */
inline DEVICE_QUALIFIER auto abs(float x) { return fabsf(x); }

/**
 * @brief Calculate the function @f$ \mathrm{sinc}(x) = \sin(\pi x)/(\pi x) @f$.
 *
 * (same convention as in @cite hockney88a). In order to avoid divisions
 * by 0, arguments whose modulus is smaller than @f$ \epsilon = 0.1 @f$
 * are evaluated by an 8th order Taylor expansion.
 * Note that the difference between sinc(x) and this expansion
 * is smaller than 0.235e-12, if x is smaller than @f$ \epsilon @f$.
 * (The next term in the expansion is the 10th order contribution, i.e.
 * @f$ \pi^{10} x^{10}/39916800 \approx 0.2346 \cdot x^{12} @f$).
 */
template <typename T> DEVICE_QUALIFIER auto sinc(T x) {
  auto constexpr epsilon = T(0.1);
#if not defined(__CUDACC__)
  using std::sin;
#endif

  auto const pix = std::numbers::pi_v<T> * x;

  if (::math::abs(x) > epsilon)
    return sin(pix) / pix;

  auto constexpr factorial = [](int n) consteval {
    int acc{1}, c{1};
    while (c < n) {
      acc *= ++c;
    }
    return acc;
  };

  /* Coefficients of the Taylor expansion of sinc */
  auto constexpr c0 = T(+1) / T(factorial(1));
  auto constexpr c2 = T(-1) / T(factorial(3));
  auto constexpr c4 = T(+1) / T(factorial(5));
  auto constexpr c6 = T(-1) / T(factorial(7));
  auto constexpr c8 = T(+1) / T(factorial(9));

  auto const pix2 = pix * pix;
  return c0 + pix2 * (c2 + pix2 * (c4 + pix2 * (c6 + pix2 * c8)));
}

/**
 * @brief One of the aliasing sums used to compute k-space errors.
 * Fortunately the one which is most important (because it converges
 * most slowly, since it is not damped exponentially) can be calculated
 * analytically. The result (which depends on the order of the spline
 * interpolation) can be written as an even trigonometric polynomial.
 * The results are tabulated here (the employed formula is eq. 7-66
 * p. 233 in @cite hockney88a).
 */
template <int cao>
DEVICE_QUALIFIER auto analytic_cotangent_sum(int n, double mesh_i) {
  static_assert(cao >= 1 and cao <= 7);
#if not defined(__CUDACC__)
  using std::cos;
#endif
  auto const theta = static_cast<double>(n) * mesh_i * std::numbers::pi;
  auto const c = Utils::sqr(cos(theta));

  if constexpr (cao == 1) {
    return 1.;
  }
  if constexpr (cao == 2) {
    return (1. + c * 2.) / 3.;
  }
  if constexpr (cao == 3) {
    return (2. + c * (11. + c * 2.)) / 15.;
  }
  if constexpr (cao == 4) {
    return (17. + c * (180. + c * (114. + c * 4.))) / 315.;
  }
  if constexpr (cao == 5) {
    return (62. + c * (1072. + c * (1452. + c * (247. + c * 2.)))) / 2835.;
  }
  if constexpr (cao == 6) {
    return (1382. +
            c * (35396. + c * (83021. + c * (34096. + c * (2026. + c * 4.))))) /
           155925.;
  }
  return (21844. +
          c * (776661. +
               c * (2801040. +
                    c * (2123860. + c * (349500. + c * (8166. + c * 4.)))))) /
         6081075.;
}

inline auto get_analytic_cotangent_sum_kernel(int cao) {
  decltype(&analytic_cotangent_sum<1>) ptr = nullptr;
  if (cao == 1) {
    ptr = &analytic_cotangent_sum<1>;
  } else if (cao == 2) {
    ptr = &analytic_cotangent_sum<2>;
  } else if (cao == 3) {
    ptr = &analytic_cotangent_sum<3>;
  } else if (cao == 4) {
    ptr = &analytic_cotangent_sum<4>;
  } else if (cao == 5) {
    ptr = &analytic_cotangent_sum<5>;
  } else if (cao == 6) {
    ptr = &analytic_cotangent_sum<6>;
  } else if (cao == 7) {
    ptr = &analytic_cotangent_sum<7>;
  }
  if (ptr == nullptr) {
    throw std::logic_error("Invalid value cao=" + std::to_string(cao));
  }
  return ptr;
}

} // namespace math
