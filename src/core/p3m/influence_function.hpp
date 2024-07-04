/*
 * Copyright (C) 2019-2024 The ESPResSo project
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

#include "p3m/common.hpp"
#include "p3m/for_each_3d.hpp"
#include "p3m/math.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <utility>
#include <vector>

/**
 * @brief Hockney/Eastwood/Ballenegger optimal influence function.
 *
 * This implements Eq. 30 of @cite cerda08d, which can be used
 * for monopole and dipole P3M by choosing the appropriate S factor.
 *
 * @tparam S Order of the differential operator, e.g. 0 for potential,
 *          1 for electric field, ...
 * @tparam m Number of aliasing terms to take into account.
 *
 * @param cao Charge assignment order.
 * @param alpha Ewald splitting parameter.
 * @param k k-vector to evaluate the function for.
 * @param h Grid spacing.
 */
template <std::size_t S, std::size_t m>
double G_opt(int cao, double alpha, Utils::Vector3d const &k,
             Utils::Vector3d const &h) {

  auto const k2 = k.norm2();
  if (k2 == 0.) {
    return 0.;
  }

  auto constexpr limit = 30.;
  auto constexpr m_start = Utils::Vector3i::broadcast(-m);
  auto constexpr m_stop = Utils::Vector3i::broadcast(m + 1);
  auto const exponent_prefactor = Utils::sqr(1. / (2. * alpha));
  auto const wavevector = (2. * std::numbers::pi) / h;
  auto const wavevector_i = 1. / wavevector;
  auto indices = Utils::Vector3i{};
  auto km = Utils::Vector3d{};
  auto fnm = Utils::Vector3d{};
  auto numerator = 0.;
  auto denominator = 0.;

  for_each_3d(
      m_start, m_stop, indices,
      [&]() {
        auto const U2 = std::pow(Utils::product(fnm), 2 * cao);
        auto const km2 = km.norm2();
        auto const exponent = exponent_prefactor * km2;
        if (exponent < limit) {
          auto const f3 = std::exp(-exponent) * (4. * std::numbers::pi / km2);
          numerator += U2 * f3 * Utils::int_pow<S>(k * km);
        }
        denominator += U2;
      },
      [&](unsigned dim, int n) {
        km[dim] = k[dim] + n * wavevector[dim];
        fnm[dim] = math::sinc(km[dim] * wavevector_i[dim]);
      });

  return numerator / (Utils::int_pow<S>(k2) * Utils::sqr(denominator));
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
 * @param params P3M parameters.
 * @param n_start Lower left corner of the grid.
 * @param n_stop Upper right corner of the grid.
 * @param KX k-space x-axis index.
 * @param KY k-space y-axis index.
 * @param KZ k-space z-axis index.
 * @param inv_box_l Inverse box length.
 * @return Values of G_opt at regular grid points.
 */
template <std::size_t S, std::size_t m = 0>
std::vector<double> grid_influence_function(P3MParameters const &params,
                                            Utils::Vector3i const &n_start,
                                            Utils::Vector3i const &n_stop,
                                            int const KX, int const KY,
                                            int const KZ,
                                            Utils::Vector3d const &inv_box_l) {

  auto const shifts = detail::calc_meshift(params.mesh);
  auto const size = n_stop - n_start;

  /* The influence function grid */
  auto g = std::vector<double>(Utils::product(size), 0.);

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (params.tuning) {
    return g;
  }

  auto const wavevector = (2. * std::numbers::pi) * inv_box_l;
  auto const half_mesh = params.mesh / 2;
  auto indices = Utils::Vector3i{};
  auto index = std::size_t(0u);

  for_each_3d(n_start, n_stop, indices, [&]() {
    if ((indices[KX] % half_mesh[0u] != 0) or
        (indices[KY] % half_mesh[1u] != 0) or
        (indices[KZ] % half_mesh[2u] != 0)) {
      auto const k =
          Utils::Vector3d{{shifts[0u][indices[KX]] * wavevector[0u],
                           shifts[1u][indices[KY]] * wavevector[1u],
                           shifts[2u][indices[KZ]] * wavevector[2u]}};
      g[index] = G_opt<S, m>(params.cao, params.alpha, k, params.a);
    }
    ++index;
  });

  return g;
}
