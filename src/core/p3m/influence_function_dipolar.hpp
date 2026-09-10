/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002-2010
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

#include <config/config.hpp>

#if defined(ESPRESSO_DP3M)

#include "p3m/common.hpp"
#include "p3m/for_each_3d.hpp"
#include "p3m/math.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

/**
 * @brief Calculate the aliasing sums for the optimal influence function.
 *
 * Evaluate the aliasing sums in the numerator and denominator of
 * the expression for the optimal influence function, see
 * @cite hockney88a eq (8-22), p. 275.
 *
 * The full equation is:
 * @f{align*}{
 * G_{\text{opt}}(\vec{n}, S, \text{cao}) &=
 *   \displaystyle\frac{1}{\sum_\vec{m} U^2}
 *   \sum_\vec{m}U^2
 *   \displaystyle\frac{\left(\vec{d}_{\text{op}}[\vec{n}](\vec{s}[\vec{n}] +
 *                            \vec{m}\odot\vec{N})\right)^S}
 *                     {\left(\vec{s}[\vec{n}] + \vec{m}\odot\vec{N}\right)^2}
 *   e^{-\pi^2\left(\alpha\vec{N}\odot\vec{a}\right)^{-2}\left(\vec{s}[\vec{n}]
 *      +\vec{m}\odot\vec{N}\right)^2} \\
 *
 * U &= \operatorname{det}\left[I_3\cdot\operatorname{sinc}\left(
 *         \vec{s}[\vec{n}]\oslash{\vec{N}}+\vec{m}\odot\vec{N}\right)
 *      \right]^{\text{cao}}
 * @f}
 *
 * with @f$ I_3 @f$ the 3x3 identity matrix, @f$ \vec{N} @f$ the global mesh
 * size in k-space coordinates, @f$ \vec{a} @f$ the mesh size, @f$ \vec{m} @f$
 * the Brillouin zone coordinates, @f$ \odot @f$ the Hadamard product,
 * @f$ \oslash @f$ the Hadamard division.
 *
 * @tparam S Order of the differential operator (2 for energy, 3 for forces)
 * @tparam m Number of Brillouin zones that contribute to the aliasing sums
 *
 * @param params      DP3M parameters
 * @param shift       shifting integer vector for a specific k-vector
 * @param d_op        differential operator for a specific k-vector
 * @return The optimal influence function for a specific k-vector.
 */
template <std::size_t S, std::size_t m>
double G_opt_dipolar(P3MParameters const &params, Utils::Vector3i const &shift,
                     Utils::Vector3i const &d_op) {

  auto constexpr exp_limit = 30.;
  auto constexpr m_start = Utils::Vector3i::broadcast(-m);
  auto constexpr m_stop = Utils::Vector3i::broadcast(m + 1);
  auto const cao = params.cao;
  auto const &mesh = params.mesh;
  auto const offset =
      Utils::hadamard_division(static_cast<Utils::Vector3d>(shift), mesh);
  auto const exp_prefactor = Utils::sqr(std::numbers::pi / params.alpha_L);
  auto indices = Utils::Vector3i{};
  auto nm = Utils::Vector3i{};
  auto fnm = Utils::Vector3d{};
  auto numerator = 0.;
  auto denominator = 0.;

  for_each_3d(
      m_start, m_stop, indices,
      [&]() {
        auto const U2 = std::pow(Utils::product(fnm), 2 * cao);
        auto const nm2 = nm.norm2();
        auto const exp_term = exp_prefactor * nm2;
        if (exp_term < exp_limit) {
          auto const f3 = U2 * std::exp(-exp_term) / nm2;
          numerator += f3 * Utils::int_pow<S>(static_cast<double>(d_op * nm));
        }
        denominator += U2;
      },
      [&](unsigned dim, int n) {
        nm[dim] = shift[dim] + n * mesh[dim];
        fnm[dim] = math::sinc(offset[dim] + n * mesh[dim]);
      });

  return numerator / (Utils::int_pow<S>(static_cast<double>(d_op.norm2())) *
                      Utils::int_pow<2>(denominator));
}

/**
 * @brief Map influence function over a grid.
 *
 * This evaluates the optimal influence function @ref G_opt_dipolar
 * over a regular grid of k vectors, and returns the values as a vector.
 *
 * @tparam S Order of the differential operator, e.g. 2 for energy, 3 for force
 * @tparam m Number of Brillouin zones that contribute to the aliasing sums
 *
 * @param params DP3M parameters
 * @param n_start Lower left corner of the grid in k-space.
 * @param n_stop Upper right corner of the grid in k-space.
 * @param inv_box_l Inverse box length
 * @return Values of the influence function at regular grid points.
 */
template <typename FloatType, std::size_t S, std::size_t m,
          Utils::MemoryOrder memory_order = Utils::MemoryOrder::ROW_MAJOR>
std::vector<FloatType> grid_influence_function_dipolar(
    P3MParameters const &params, Utils::Vector3i const &n_start,
    Utils::Vector3i const &n_stop, Utils::Vector3d const &inv_box_l) {

  auto const size = n_stop - n_start;

  /* The influence function grid */
  auto g = std::vector<FloatType>(Utils::product(size), FloatType(0));

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (params.tuning) {
    return g;
  }

  auto prefactor =
      Utils::product(params.mesh) * 2. * Utils::int_pow<2>(inv_box_l[0]);

  auto const offset = calc_p3m_mesh_shift(params.mesh, false)[0];
  auto const d_op = calc_p3m_mesh_shift(params.mesh, true)[0];
  auto const half_mesh = params.mesh / 2;
  auto indices = Utils::Vector3i{};
  auto shift_off = Utils::Vector3i{};
  auto d_op_off = Utils::Vector3i{};

  for_each_3d(
      n_start, n_stop, indices,
      [&]() {
        if (((indices[0u] % half_mesh[0u] != 0) or
             (indices[1u] % half_mesh[1u] != 0) or
             (indices[2u] % half_mesh[2u] != 0))) {
          auto const index =
              Utils::get_linear_index<memory_order>(indices - n_start, size);
          g[index] = FloatType(
              prefactor * G_opt_dipolar<S, m>(params, shift_off, d_op_off));
        }
      },
      [&](unsigned dim, int n) {
        d_op_off[dim] = d_op[n];
        shift_off[dim] = offset[n];
      });

  return g;
}

template <std::size_t m>
inline double G_opt_dipolar_self_energy(P3MParameters const &params,
                                        Utils::Vector3i const &shift) {

  auto constexpr m_start = Utils::Vector3i::broadcast(-m);
  auto constexpr m_stop = Utils::Vector3i::broadcast(m + 1);
  auto const cao = params.cao;
  auto const &mesh = params.mesh;
  auto const offset =
      Utils::hadamard_division(static_cast<Utils::Vector3d>(shift), mesh);
  auto indices = Utils::Vector3i{};
  auto fnm = Utils::Vector3d{};
  auto energy = 0.;

  for_each_3d(
      m_start, m_stop, indices,
      [&]() { energy += std::pow(Utils::product(fnm), 2 * cao); },
      [&](unsigned dim, int n) {
        fnm[dim] = math::sinc(offset[dim] + n * mesh[dim]);
      });

  return energy;
}

/**
 * @brief Calculate self-energy of the influence function.
 *
 * @tparam m Number of Brillouin zones that contribute to the aliasing sums
 *
 * @param params DP3M parameters
 * @param n_start Lower left corner of the grid in k-space.
 * @param n_stop Upper right corner of the grid in k-space.
 * @param g Energies on the grid.
 * @return Total self-energy.
 */
template <typename FloatType, std::size_t m>
inline double grid_influence_function_self_energy(
    P3MParameters const &params, Utils::Vector3i const &n_start,
    Utils::Vector3i const &n_stop, std::vector<FloatType> const &g) {

  auto const offset = calc_p3m_mesh_shift(params.mesh, false)[0];
  auto const d_op = calc_p3m_mesh_shift(params.mesh, true)[0];
  auto const half_mesh = params.mesh / 2;
  auto indices = Utils::Vector3i{};
  auto shift_off = Utils::Vector3i{};
  auto d_op_off = Utils::Vector3i{};
  auto index = std::size_t(0u);
  auto energy = 0.;

  for_each_3d(
      n_start, n_stop, indices,
      [&]() {
        if (((indices[0u] % half_mesh[0u] != 0) or
             (indices[1u] % half_mesh[1u] != 0) or
             (indices[2u] % half_mesh[2u] != 0))) {
          auto const U2 = G_opt_dipolar_self_energy<m>(params, shift_off);
          energy += static_cast<double>(g[index]) * U2 * d_op_off.norm2();
        }
        ++index;
      },
      [&](unsigned dim, int n) {
        d_op_off[dim] = d_op[n];
        shift_off[dim] = offset[n];
      });

  return energy;
}

#endif // defined(ESPRESSO_DP3M)
