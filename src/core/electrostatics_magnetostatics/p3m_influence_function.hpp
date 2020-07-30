/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#ifndef ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP
#define ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>

/**
 * @brief Ewald method Green's function in k space.
 *
 * Really this is just the potential of a point charge
 * at the origin screened by Gaussian charge density of
 * opposite sign. The Ewald parameter is really just the
 * width of the Gaussian.
 *
 * @tparam T Floating-point type
 * @param alpha Ewald parameter.
 * @param k2 Square of magnitude of the k vector.
 * @return Green's function evaluated at k2
 */
template <typename T> T g_ewald(T alpha, T k2) {
  auto constexpr limit = T{30};
  auto const exponent = Utils::sqr(1. / (2. * alpha)) * k2;
  return (exponent < limit) ? std::exp(-exponent) * (4. * Utils::pi<T>() / k2)
                            : 0.0;
}

template <size_t cao, unsigned S, unsigned m = 0>
inline std::pair<double, double> aliasing_sums_ik(double alpha,
                                                  const Utils::Vector3d &k,
                                                  const Utils::Vector3d &h) {
  using Utils::int_pow;
  using Utils::pi;
  using Utils::sinc;
  using Utils::Vector3d;

  constexpr int RX = 0;
  constexpr int RY = 1;
  constexpr int RZ = 2;

  double numerator = 0.0;
  double denominator = 0.0;

  for (int mx = -m; mx <= m; mx++) {
    for (int my = -m; my <= m; my++) {
      for (int mz = -m; mz <= m; mz++) {
        auto const km =
            k + 2 * pi() * Vector3d{mx / h[RX], my / h[RY], mz / h[RZ]};
        auto const U2 = int_pow<2 * cao>(sinc(0.5 * km[RX] * h[RX] / pi()) *
                                         sinc(0.5 * km[RY] * h[RY] / pi()) *
                                         sinc(0.5 * km[RZ] * h[RZ] / pi()));

        numerator += U2 * g_ewald(alpha, km.norm2()) * int_pow<S>(k * km);
        denominator += U2;
      }
    }
  }
  return {numerator, denominator};
}

template <size_t cao, size_t S, size_t m = 0>
double G_opt(double alpha, const Utils::Vector3d &k, const Utils::Vector3d &h) {
  using Utils::int_pow;
  using Utils::sqr;

  auto const k2 = k.norm2();
  if (k2 == 0.0) {
    return 0.0;
  }

  const auto as = aliasing_sums_ik<cao, S, m>(alpha, k, h);
  return as.first / (int_pow<S>(k2) * sqr(as.second));
}

#endif // ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP
