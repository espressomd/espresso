/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef ESPRESSO_DP3M_INFLUENCE_FUNCTION_HPP
#define ESPRESSO_DP3M_INFLUENCE_FUNCTION_HPP

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
#include <vector>

#if defined(DP3M)

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
  using Utils::int_pow;
  using Utils::sinc;
  constexpr double limit = 30;

  double numerator = 0.0;
  double denominator = 0.0;

  auto const f1 = 1.0 / static_cast<double>(params.mesh[0]);
  auto const f2 = Utils::sqr(Utils::pi() / params.alpha_L);

  for (double mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    auto const nmx = shift[0] + params.mesh[0] * mx;
    auto const sx = std::pow(sinc(f1 * nmx), 2.0 * params.cao);
    for (double my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      auto const nmy = shift[1] + params.mesh[0] * my;
      auto const sy = sx * std::pow(sinc(f1 * nmy), 2.0 * params.cao);
      for (double mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        auto const nmz = shift[2] + params.mesh[0] * mz;
        auto const sz = sy * std::pow(sinc(f1 * nmz), 2.0 * params.cao);
        auto const nm2 = Utils::sqr(nmx) + Utils::sqr(nmy) + Utils::sqr(nmz);
        auto const exponent = f2 * nm2;
        if (exponent < limit) {
          auto const f3 = sz * std::exp(-exponent) / nm2;
          auto const n_nm = d_op[0] * nmx + d_op[1] * nmy + d_op[2] * nmz;
          numerator += f3 * int_pow<S>(n_nm);
        }
        denominator += sz;
      }
    }
  }
  return numerator / (int_pow<S>(static_cast<double>(d_op.norm2())) *
                      Utils::sqr(denominator));
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
 * @param n_end Upper right corner of the grid.
 * @param box_l Box size
 * @return Values of the influence function at regular grid points.
 */
template <std::size_t S>
std::vector<double> grid_influence_function(P3MParameters const &params,
                                            Utils::Vector3i const &n_start,
                                            Utils::Vector3i const &n_end,
                                            Utils::Vector3d const &box_l) {

  auto const size = n_end - n_start;

  /* The influence function grid */
  auto g =
      std::vector<double>(boost::accumulate(size, 1, std::multiplies<>()), 0.);

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (params.tuning) {
    return g;
  }

  double fak1 = Utils::int_pow<3>(static_cast<double>(params.mesh[0])) * 2.0 /
                Utils::sqr(box_l[0]);

  auto const shifts = detail::calc_meshift(params.mesh, false);
  auto const d_ops = detail::calc_meshift(params.mesh, true);

  Utils::Vector3i n{};
  for (n[0] = n_start[0]; n[0] < n_end[0]; n[0]++) {
    for (n[1] = n_start[1]; n[1] < n_end[1]; n[1]++) {
      for (n[2] = n_start[2]; n[2] < n_end[2]; n[2]++) {
        auto const ind = Utils::get_linear_index(n - n_start, size,
                                                 Utils::MemoryOrder::ROW_MAJOR);

        if (((n[0] % (params.mesh[0] / 2) == 0) &&
             (n[1] % (params.mesh[0] / 2) == 0) &&
             (n[2] % (params.mesh[0] / 2) == 0))) {
          g[ind] = 0.0;
        } else {
          auto const shift = Utils::Vector3i{shifts[0][n[0]], shifts[0][n[1]],
                                             shifts[0][n[2]]};
          auto const d_op =
              Utils::Vector3i{d_ops[0][n[0]], d_ops[0][n[1]], d_ops[0][n[2]]};
          auto const fak2 = G_opt_dipolar<S>(params, shift, d_op);
          g[ind] = fak1 * fak2;
        }
      }
    }
  }
  return g;
}

double G_opt_dipolar_self_energy(P3MParameters const &params,
                                 Utils::Vector3i const &shift) {
  using Utils::sinc;
  double u_sum = 0.0;
  constexpr int limit = P3M_BRILLOUIN + 5;

  auto const f1 = 1.0 / static_cast<double>(params.mesh[0]);

  for (double mx = -limit; mx <= limit; mx++) {
    auto const nmx = shift[0] + params.mesh[0] * mx;
    auto const sx = std::pow(sinc(f1 * nmx), 2.0 * params.cao);
    for (double my = -limit; my <= limit; my++) {
      auto const nmy = shift[1] + params.mesh[0] * my;
      auto const sy = sx * std::pow(sinc(f1 * nmy), 2.0 * params.cao);
      for (double mz = -limit; mz <= limit; mz++) {
        auto const nmz = shift[2] + params.mesh[0] * mz;
        auto const sz = sy * std::pow(sinc(f1 * nmz), 2.0 * params.cao);
        u_sum += sz;
      }
    }
  }
  return u_sum;
}

/**
 * @brief Calculate self-energy of the influence function.
 *
 * @param params DP3M parameters
 * @param n_start Lower left corner of the grid
 * @param n_end Upper right corner of the grid.
 * @param g Energies on the grid.
 * @return Total self-energy.
 */
double grid_influence_function_self_energy(P3MParameters const &params,
                                           Utils::Vector3i const &n_start,
                                           Utils::Vector3i const &n_end,
                                           std::vector<double> const &g) {
  auto const size = n_end - n_start;

  auto const shifts = detail::calc_meshift(params.mesh, false);
  auto const d_ops = detail::calc_meshift(params.mesh, true);

  double energy = 0.0;
  Utils::Vector3i n{};
  for (n[0] = n_start[0]; n[0] < n_end[0]; n[0]++) {
    for (n[1] = n_start[1]; n[1] < n_end[1]; n[1]++) {
      for (n[2] = n_start[2]; n[2] < n_end[2]; n[2]++) {
        if (((n[0] % (params.mesh[0] / 2) == 0) &&
             (n[1] % (params.mesh[0] / 2) == 0) &&
             (n[2] % (params.mesh[0] / 2) == 0))) {
          energy += 0.0;
        } else {
          auto const ind = Utils::get_linear_index(
              n - n_start, size, Utils::MemoryOrder::ROW_MAJOR);
          auto const shift = Utils::Vector3i{shifts[0][n[0]], shifts[0][n[1]],
                                             shifts[0][n[2]]};
          auto const d_op =
              Utils::Vector3i{d_ops[0][n[0]], d_ops[0][n[1]], d_ops[0][n[2]]};
          auto const U2 = G_opt_dipolar_self_energy(params, shift);
          energy += g[ind] * U2 * d_op.norm2();
        }
      }
    }
  }
  return energy;
}

#endif
#endif
