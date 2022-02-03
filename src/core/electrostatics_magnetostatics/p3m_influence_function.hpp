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

#include "electrostatics_magnetostatics/p3m-common.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/index.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>

#include <boost/range/numeric.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace detail {
template <typename T> T g_ewald(T alpha, T k2) {
  auto constexpr limit = T{30};
  auto const exponent = Utils::sqr(1. / (2. * alpha)) * k2;
  return (exponent < limit) ? std::exp(-exponent) * (4. * Utils::pi<T>() / k2)
                            : 0.0;
}

template <std::size_t S, std::size_t m>
std::pair<double, double> aliasing_sums_ik(std::size_t cao, double alpha,
                                           const Utils::Vector3d &k,
                                           const Utils::Vector3d &h) {
  using namespace detail::FFT_indexing;
  using Utils::int_pow;
  using Utils::sinc;
  using Utils::Vector3d;

  constexpr double two_pi = 2 * Utils::pi();
  constexpr double two_pi_i = 1 / two_pi;

  double numerator = 0.0;
  double denominator = 0.0;

  for (int mx = -m; mx <= m; mx++) {
    for (int my = -m; my <= m; my++) {
      for (int mz = -m; mz <= m; mz++) {
        auto const km =
            k + two_pi * Vector3d{mx / h[RX], my / h[RY], mz / h[RZ]};
        auto const U2 = std::pow(sinc(km[RX] * h[RX] * two_pi_i) *
                                     sinc(km[RY] * h[RY] * two_pi_i) *
                                     sinc(km[RZ] * h[RZ] * two_pi_i),
                                 2 * cao);

        numerator += U2 * g_ewald(alpha, km.norm2()) * int_pow<S>(k * km);
        denominator += U2;
      }
    }
  }
  return {numerator, denominator};
}
} // namespace detail

/**
 * @brief Hockney/Eastwood/Ballenegger optimal influence function.
 *
 *  This implements Eq. 30 of @cite cerda08d, which can be used
 *  for monopole and dipole P3M by choosing the appropriate S factor.
 *
 * @tparam S Order of the differential operator, e.g. 0 for potential,
 *          1 for electric field, ...
 * @tparam m Number of aliasing terms to take into account.
 * @tparam T Floating-point type.
 *
 * @param cao Charge assignment order.
 * @param alpha Ewald splitting parameter.
 * @param k k Vector to evaluate the function for.
 * @param h Grid spacing.
 */
template <std::size_t S, std::size_t m, class T>
double G_opt(std::size_t cao, T alpha, const Utils::Vector3<T> &k,
             const Utils::Vector3<T> &h) {
  using Utils::int_pow;
  using Utils::sqr;

  auto const k2 = k.norm2();
  if (k2 == 0.0) {
    return 0.0;
  }

  const auto as = detail::aliasing_sums_ik<S, m>(cao, alpha, k, h);
  return as.first / (int_pow<S>(k2) * sqr(as.second));
}

/**
 * @brief Map influence function over a grid.
 *
 * This evaluates the optimal influence function @ref G_opt
 * over a regular grid of k vectors, and returns the values as a vector.
 *
 * @tparam S Order of the differential operator, e.g. 0 for potential,
 *          1 for electric field...
 * @tparam m Number of aliasing terms to take into account.
 *
 * @param params P3M parameters
 * @param n_start Lower left corner of the grid
 * @param n_end Upper right corner of the grid.
 * @param box_l Box size
 * @return Values of G_opt at regular grid points.
 */
template <std::size_t S, std::size_t m = 0>
std::vector<double> grid_influence_function(const P3MParameters &params,
                                            const Utils::Vector3i &n_start,
                                            const Utils::Vector3i &n_end,
                                            const Utils::Vector3d &box_l) {
  using namespace detail::FFT_indexing;

  auto const shifts = detail::calc_meshift(params.mesh);

  auto const size = n_end - n_start;

  /* The influence function grid */
  auto g =
      std::vector<double>(boost::accumulate(size, 1, std::multiplies<>()), 0.);

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (params.tuning) {
    return g;
  }

  auto const h = Utils::Vector3d{params.a};

  Utils::Vector3i n{};
  for (n[0] = n_start[0]; n[0] < n_end[0]; n[0]++) {
    for (n[1] = n_start[1]; n[1] < n_end[1]; n[1]++) {
      for (n[2] = n_start[2]; n[2] < n_end[2]; n[2]++) {
        auto const ind = Utils::get_linear_index(n - n_start, size,
                                                 Utils::MemoryOrder::ROW_MAJOR);
        if ((n[KX] % (params.mesh[RX] / 2) == 0) &&
            (n[KY] % (params.mesh[RY] / 2) == 0) &&
            (n[KZ] % (params.mesh[RZ] / 2) == 0)) {
          g[ind] = 0.0;
        } else {
          auto const k = 2 * Utils::pi() *
                         Utils::Vector3d{shifts[RX][n[KX]] / box_l[RX],
                                         shifts[RY][n[KY]] / box_l[RY],
                                         shifts[RZ][n[KZ]] / box_l[RZ]};

          g[ind] = G_opt<S, m>(params.cao, params.alpha, k, h);
        }
      }
    }
  }

  return g;
}

#endif // ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP
