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

#include <config/config.hpp>

#if defined(ESPRESSO_P3M)

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
 * @brief Calculate the aliasing sums for the optimal influence function.
 *
 * This implements Eq. 30 of @cite cerda08d, which can be used
 * for monopole and dipole P3M by choosing the appropriate S factor.
 * @f[
 * \tilde{G}_{\text{opt}} ( k ) = \frac{\sum_{m \in \mathbb{Z}^{3}} [ [ \tilde{
 * \left ( D \right ) } ( k ) \cdot i k_{m} ]^{S} ( \hat{U} ( k_{m} ) )^{2}
 * \hat{\phi} ( k_{m} ) ]}{[ \tilde{ \left ( D \right ) } ( k ) ]^{2 S} [
 * \sum_{m \in \mathbb{Z}^{3}} ( \hat{U} ( k_{m} ) )^{2} ]^{2}}
 * @f]
 *
 * Eq 8.29 in Hockney:
 * @f[
 * G_{\text{opt}}(\mathbf{k}) = \frac{\hat{\mathbf{D}}\sum_n
 * \hat{\mathbf{R}}^\ast \hat{U}^2}{|\hat{\mathbf{D}}|^2\sum_n
 * \hat{U}^2\sum_{n'}^{'} \hat{U}^2}
 * @f]
 *
 * The full equation is:
 * @f{align*}{
 * G_{\text{opt}}(\vec{n}, S, \text{cao}) &=
 *   \displaystyle\frac{4\pi}{\sum_\vec{m} U^2}
 *   \sum_\vec{m}U^2
 *   \displaystyle\frac{\left(\vec{k}\odot\vec{k}_m\right)^S}
 *                     {|\vec{k}_m|^2}
 *   e^{-1/(2\alpha)^2|\vec{k}_m|^2} \\
 *
 * U &= \operatorname{det}\left[I_3\cdot\operatorname{sinc}\left(
 *          \frac{\vec{k}_m\odot\vec{a}}{2\pi}\right)\right]^{\text{cao}} \\
 *
 * \vec{k} &= \frac{2\pi}{\vec{N}\odot\vec{a}}\vec{s}[\vec{n}] \\
 *
 * \vec{k}_m &= \frac{2\pi}{\vec{N}\odot\vec{a}}
 *              \left(\vec{s}[\vec{n}]+\vec{m}\odot\vec{N}\right)
 * @f}
 *
 * with @f$ I_3 @f$ the 3x3 identity matrix, @f$ \vec{N} @f$ the global mesh
 * size in k-space coordinates, @f$ \vec{a} @f$ the mesh size, @f$ \vec{m} @f$
 * the Brillouin zone coordinates, @f$ \odot @f$ the Hadamard product.
 *
 * @tparam S Order of the differential operator (0 for potential, 1 for E-field)
 * @tparam m Number of Brillouin zones that contribute to the aliasing sums
 *
 * @param params   P3M parameters
 * @param k        k-vector to evaluate the function for.
 * @return The optimal influence function for a specific k-vector.
 */
template <std::size_t S, std::size_t m>
double G_opt(P3MParameters const &params, Utils::Vector3d const &k) {

  auto const k2 = k.norm2();
  if (k2 == 0.) {
    return 0.;
  }

  auto constexpr exp_limit = 30.;
  auto constexpr m_start = Utils::Vector3i::broadcast(-m);
  auto constexpr m_stop = Utils::Vector3i::broadcast(m + 1);
  auto const cao = params.cao;
  auto const exponent_prefactor = Utils::sqr(1. / (2. * params.alpha));
  auto const wavevector = (2. * std::numbers::pi) / params.a;
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
        auto const exp_term = exponent_prefactor * km2;
        if (exp_term < exp_limit) {
          auto const f3 = std::exp(-exp_term) * (4. * std::numbers::pi / km2);
          numerator += U2 * f3 * Utils::int_pow<S>(k * km);
        }
        denominator += U2;
      },
      [&](unsigned dim, int n) {
        km[dim] = k[dim] + n * wavevector[dim];
        fnm[dim] = math::sinc(km[dim] * wavevector_i[dim]);
      });

  if (numerator == 0. and denominator != 0.) {
    // avoid FPE issues involving subnormal numbers
    return 0.;
  }
  return numerator / (Utils::int_pow<S>(k2) * Utils::sqr(denominator));
}

/**
 * @brief Map influence function over a grid.
 *
 * This evaluates the optimal influence function @ref G_opt
 * over a regular grid of k vectors, and returns the values as a vector.
 *
 * @tparam S Order of the differential operator (0 for potential, 1 for E-field)
 * @tparam m Number of Brillouin zones that contribute to the aliasing sums
 *
 * @param params P3M parameters.
 * @param n_start Lower left corner of the grid in k-space.
 * @param n_stop Upper right corner of the grid in k-space.
 * @param inv_box_l Inverse box length.
 * @return Values of G_opt at regular grid points.
 */
template <typename FloatType, std::size_t S, std::size_t m,
          Utils::MemoryOrder memory_order>
std::vector<FloatType> grid_influence_function(
    P3MParameters const &params, Utils::Vector3i const &n_start,
    Utils::Vector3i const &n_stop, Utils::Vector3d const &inv_box_l) {

  auto const shifts = calc_p3m_mesh_shift(params.mesh);
  auto const size = n_stop - n_start;

  /* The influence function grid */
  auto g = std::vector<FloatType>(Utils::product(size), FloatType(0));

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (params.tuning) {
    return g;
  }

  auto const wavevector = (2. * std::numbers::pi) * inv_box_l;
  auto const half_mesh = params.mesh / 2;
  auto indices = Utils::Vector3i{};

  for_each_3d(n_start, n_stop, indices, [&]() {
    if ((indices[0u] % half_mesh[0u] != 0) or
        (indices[1u] % half_mesh[1u] != 0) or
        (indices[2u] % half_mesh[2u] != 0)) {
      auto const k =
          Utils::Vector3d{{shifts[0u][indices[0u]] * wavevector[0u],
                           shifts[1u][indices[1u]] * wavevector[1u],
                           shifts[2u][indices[2u]] * wavevector[2u]}};
      auto const index =
          Utils::get_linear_index<memory_order>(indices - n_start, size);
      g[index] = FloatType(G_opt<S, m>(params, k));
    }
  });

  return g;
}

#endif // defined(ESPRESSO_P3M)
