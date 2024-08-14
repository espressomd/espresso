/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#if defined(DP3M)

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

/** Calculate the aliasing sums for the optimal influence function.
 *
 *  Calculates the aliasing sums in the numerator and denominator of
 *  the expression for the optimal influence function (see
 *  @cite hockney88a : 8-22, p. 275).
 *
 *  \tparam S          order (2 for energy, 3 for forces)
 *  \param params      DP3M parameters
 *  \param shift       shift for a given n-vector
 *  \param d_op        differential operator for a given n-vector
 *  \return The result of the fraction.
 */
template <std::size_t S>
double G_opt_dipolar(P3MParameters const &params, Utils::Vector3i const &shift,
                     Utils::Vector3i const &d_op) {

  auto constexpr limit = P3M_BRILLOUIN;
  auto constexpr exp_limit = 30.;
  auto constexpr m_start = Utils::Vector3i::broadcast(-limit);
  auto constexpr m_stop = Utils::Vector3i::broadcast(limit + 1);
  auto const cao = params.cao;
  auto const mesh = params.mesh[0];
  auto const offset =
      static_cast<Utils::Vector3d>(shift) / static_cast<double>(mesh);
  auto const f2 = Utils::sqr(std::numbers::pi / params.alpha_L);
  auto indices = Utils::Vector3i{};
  auto nm = Utils::Vector3i{};
  auto fnm = Utils::Vector3d{};
  auto numerator = 0.;
  auto denominator = 0.;

  for_each_3d(
      m_start, m_stop, indices,
      [&]() {
        auto const norm_sq = nm.norm2();
        auto const sz = std::pow(Utils::product(fnm), 2 * cao);
        auto const exp_term = f2 * norm_sq;
        if (exp_term < exp_limit) {
          auto const f3 = sz * std::exp(-exp_term) / norm_sq;
          numerator += f3 * Utils::int_pow<S>(d_op * nm);
        }
        denominator += sz;
      },
      [&](unsigned dim, int n) {
        nm[dim] = shift[dim] + n * mesh;
        fnm[dim] = math::sinc(offset[dim] + n * mesh);
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
 *
 * @param params DP3M parameters
 * @param n_start Lower left corner of the grid
 * @param n_stop Upper right corner of the grid.
 * @param inv_box_l Inverse box length
 * @return Values of the influence function at regular grid points.
 */
template <typename FloatType, std::size_t S>
std::vector<FloatType> grid_influence_function(
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

  auto prefactor = Utils::int_pow<3>(static_cast<double>(params.mesh[0])) * 2. *
                   Utils::int_pow<2>(inv_box_l[0]);

  auto const offset = detail::calc_meshift(params.mesh, false)[0];
  auto const d_op = detail::calc_meshift(params.mesh, true)[0];
  auto const half_mesh = params.mesh[0] / 2;
  auto indices = Utils::Vector3i{};
  auto shift_off = Utils::Vector3i{};
  auto d_op_off = Utils::Vector3i{};
  auto index = std::size_t(0u);

  for_each_3d(
      n_start, n_stop, indices,
      [&]() {
        if (((indices[0] % half_mesh != 0) or (indices[1] % half_mesh != 0) or
             (indices[2] % half_mesh != 0))) {
          g[index] = FloatType(prefactor *
                               G_opt_dipolar<S>(params, shift_off, d_op_off));
        }
        ++index;
      },
      [&](unsigned dim, int n) {
        d_op_off[dim] = d_op[n];
        shift_off[dim] = offset[n];
      });

  return g;
}

inline double G_opt_dipolar_self_energy(P3MParameters const &params,
                                        Utils::Vector3i const &shift) {

  auto constexpr limit = P3M_BRILLOUIN + 1;
  auto constexpr m_start = Utils::Vector3i::broadcast(-limit);
  auto constexpr m_stop = Utils::Vector3i::broadcast(limit + 1);
  auto const cao = params.cao;
  auto const mesh = params.mesh[0];
  auto const offset =
      static_cast<Utils::Vector3d>(shift) / static_cast<double>(mesh);
  auto indices = Utils::Vector3i{};
  auto fnm = Utils::Vector3d{};
  auto energy = 0.;

  for_each_3d(
      m_start, m_stop, indices,
      [&]() { energy += std::pow(Utils::product(fnm), 2 * cao); },
      [&](unsigned dim, int n) {
        fnm[dim] = math::sinc(offset[dim] + n * mesh);
      });

  return energy;
}

/**
 * @brief Calculate self-energy of the influence function.
 *
 * @param params DP3M parameters
 * @param n_start Lower left corner of the grid
 * @param n_stop Upper right corner of the grid.
 * @param g Energies on the grid.
 * @return Total self-energy.
 */
template <typename FloatType>
inline double grid_influence_function_self_energy(
    P3MParameters const &params, Utils::Vector3i const &n_start,
    Utils::Vector3i const &n_stop, std::vector<FloatType> const &g) {

  auto const offset = detail::calc_meshift(params.mesh, false)[0];
  auto const d_op = detail::calc_meshift(params.mesh, true)[0];
  auto const half_mesh = params.mesh[0] / 2;
  auto indices = Utils::Vector3i{};
  auto shift_off = Utils::Vector3i{};
  auto d_op_off = Utils::Vector3i{};
  auto index = std::size_t(0u);
  auto energy = 0.;

  for_each_3d(
      n_start, n_stop, indices,
      [&]() {
        if (((indices[0] % half_mesh != 0) or (indices[1] % half_mesh != 0) or
             (indices[2] % half_mesh != 0))) {
          auto const U2 = G_opt_dipolar_self_energy(params, shift_off);
          energy += double(g[index]) * U2 * d_op_off.norm2();
        }
        ++index;
      },
      [&](unsigned dim, int n) {
        d_op_off[dim] = d_op[n];
        shift_off[dim] = offset[n];
      });

  return energy;
}

#endif // defined(DP3M)
