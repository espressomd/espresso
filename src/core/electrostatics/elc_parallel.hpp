/*
 * Copyright (C) 2026 The ESPResSo project
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

/** @file
 *  Shared-memory parallel helpers for the ELC kernels.
 */

#include <Kokkos_Core.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace ELC {

/** @brief Kokkos reducer summing a fixed-size array of doubles. */
template <std::size_t N, class Kernel> struct ArraySum {
  using value_type = std::array<double, N>;
  Kernel kernel;

  void operator()(std::size_t i, value_type &acc) const { kernel(i, acc); }
  void init(value_type &acc) const { acc.fill(0.); }
  void join(value_type &dst, value_type const &src) const {
    for (std::size_t k = 0; k < N; ++k) {
      dst[k] += src[k];
    }
  }
};

/**
 * @brief Sum per-particle contributions into a fixed-size array.
 * @param name    Kokkos kernel label
 * @param n       number of particles
 * @param kernel  callable ``(std::size_t i, std::array<double, N> &acc)``
 *                adding the contribution of particle @p i to @p acc
 */
template <std::size_t N>
std::array<double, N> parallel_array_sum(char const *name, std::size_t n,
                                         auto const &kernel) {
  std::array<double, N> result{};
  if (Kokkos::num_threads() == 1) {
    for (std::size_t i = 0; i < n; ++i) {
      kernel(i, result);
    }
    return result;
  }
  using Reducer = ArraySum<N, std::remove_cvref_t<decltype(kernel)>>;
  Kokkos::parallel_reduce(
      name, Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
      Reducer{kernel}, result);
  return result;
}

/**
 * @brief One accumulation row per host thread, folded in thread order.
 * The fold order is fixed, so the result is deterministic for a given
 * number of threads.
 */
class PerThreadRows {
  std::size_t m_n_cols;
  std::size_t m_stride;
  Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace> m_rows;

public:
  explicit PerThreadRows(std::size_t n_cols)
      : m_n_cols{n_cols}, m_stride{(n_cols + 7u) / 8u * 8u},
        m_rows{"ELC::PerThreadRows",
               static_cast<std::size_t>(
                   Kokkos::DefaultHostExecutionSpace().concurrency()),
               m_stride} {}

  double *row(int tid) {
    assert(static_cast<std::size_t>(tid) < m_rows.extent(0));
    return &m_rows(static_cast<std::size_t>(tid), 0);
  }

  std::vector<double> fold() const {
    std::vector<double> result(m_n_cols, 0.);
    for (std::size_t tid = 0; tid < m_rows.extent(0); ++tid) {
      for (std::size_t k = 0; k < m_n_cols; ++k) {
        result[k] += m_rows(tid, k);
      }
    }
    return result;
  }
};

} // namespace ELC
