/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#include <utils/index.hpp>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#endif

#include <concepts>
#include <cstddef>
#include <type_traits>

namespace detail {

constexpr inline void noop_projector(unsigned, int) {}

template <typename T>
concept IndexVectorConcept = requires(T vector) {
  { vector[0] } -> std::convertible_to<std::size_t>;
};

} // namespace detail

/**
 * @brief Repeat an operation on every element of a 3D grid.
 *
 * Intermediate values that depend on the iterated coordinates
 * are calculated and stored once per iteration. This is useful
 * when the operation is costly.
 *
 * @param start       Initial values for the loop counters.
 * @param stop        Final values (one-past-the-end) for the loop counters.
 * @param counters    Loop counters.
 * @param kernel      Functor to execute.
 * @param projector   Projection of the current loop counter.
 * @tparam Kernel     Nullary function.
 * @tparam Projector  Binary function that takes a nesting depth and a loop
 *                    counter as arguments and projects a value.
 */
template <class Kernel, class Projector = decltype(detail::noop_projector)>
  requires std::invocable<Kernel> and std::invocable<Projector, unsigned, int>
void for_each_3d(detail::IndexVectorConcept auto &&start,
                 detail::IndexVectorConcept auto &&stop,
                 detail::IndexVectorConcept auto &&counters, Kernel &&kernel,
                 Projector &&projector = detail::noop_projector) {
  auto &nx = counters[0u];
  auto &ny = counters[1u];
  auto &nz = counters[2u];
  for (nx = start[0u]; nx < stop[0u]; ++nx) {
    projector(0u, nx);
    for (ny = start[1u]; ny < stop[1u]; ++ny) {
      projector(1u, ny);
      for (nz = start[2u]; nz < stop[2u]; ++nz) {
        projector(2u, nz);
        kernel();
      }
    }
  }
}

/**
 * @brief Repeat an operation on every element of a 3D grid.
 *
 * Intermediate values that depend on the iterated coordinates
 * are calculated and stored once per iteration. This is useful
 * when the operation is costly.
 *
 * @param start       Initial values for the loop counters.
 * @param stop        Final values (one-past-the-end) for the loop counters.
 * @param counters    Loop counters.
 * @param kernel      Functor to execute.
 * @param projector   Projection of the current loop counter.
 * @tparam Order      Data layout.
 * @tparam Kernel     Nullary function.
 * @tparam Projector  Binary function that takes a nesting depth and a loop
 *                    counter as arguments and projects a value.
 */
template <Utils::MemoryOrder Order, class Kernel,
          class Projector = decltype(detail::noop_projector)>
  requires std::invocable<Kernel> and std::invocable<Projector, unsigned, int>
void for_each_3d_order(detail::IndexVectorConcept auto &&start,
                       detail::IndexVectorConcept auto &&stop,
                       detail::IndexVectorConcept auto &&counters,
                       Kernel &&kernel,
                       Projector &&projector = detail::noop_projector) {
  auto constexpr is_row_major = Order == Utils::MemoryOrder::ROW_MAJOR;
  auto constexpr index_fast = is_row_major ? 2u : 0u;
  auto constexpr index_slow = is_row_major ? 0u : 2u;
  auto constexpr index_medium = 1u;
  auto &nx = counters[index_slow];
  auto &ny = counters[index_medium];
  auto &nz = counters[index_fast];
  for (nx = start[index_slow]; nx < stop[index_slow]; ++nx) {
    projector(index_slow, nx);
    for (ny = start[index_medium]; ny < stop[index_medium]; ++ny) {
      projector(index_medium, ny);
      for (nz = start[index_fast]; nz < stop[index_fast]; ++nz) {
        projector(index_fast, nz);
        kernel();
      }
    }
  }
}

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
/** @brief Mapping between ESPResSo and Kokkos tags for memory order */
template <Utils::MemoryOrder Order>
using LayoutIterate = std::conditional_t<
    Order == Utils::MemoryOrder::COLUMN_MAJOR,
    std::integral_constant<Kokkos::Iterate, Kokkos::Iterate::Left>,
    std::integral_constant<Kokkos::Iterate, Kokkos::Iterate::Right>>;
#endif

/**
 * @brief Run a kernel(index_3d, linear_index) over the given 3d range with
 * given memory order
 */
template <Utils::MemoryOrder memory_order, class Kernel>
void for_each_3d_lin(detail::IndexVectorConcept auto &&start,
                     detail::IndexVectorConcept auto &&stop, Kernel &&kernel) {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  if (Kokkos::num_threads() > 1) {
    auto const size = stop - start;
    constexpr Kokkos::Iterate iter = LayoutIterate<memory_order>::value;
    using Range3d = Kokkos::MDRangePolicy<Kokkos::Rank<3, iter, iter>>;
    Range3d policy({0, 0, 0}, {size[0], size[1], size[2]});
    Kokkos::parallel_for(
        "for_each_3d", policy, KOKKOS_LAMBDA(int i, int j, int k) {
          auto const linear_idx =
              Utils::get_linear_index<memory_order>(i, j, k, size);
          kernel(start + Utils::Vector3i{{i, j, k}}, linear_idx);
        });
    return;
  }
#endif

  int linear_loop_index = 0u;
  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    for (int nx = start[0u]; nx < stop[0u]; ++nx) {
      for (int ny = start[1u]; ny < stop[1u]; ++ny) {
        for (int nz = start[2u]; nz < stop[2u]; ++nz) {
          kernel(Utils::Vector3i{nx, ny, nz}, linear_loop_index);
          linear_loop_index++;
        }
      }
    }
  } else {
    for (int nz = start[2u]; nz < stop[2u]; ++nz) {
      for (int ny = start[1u]; ny < stop[1u]; ++ny) {
        for (int nx = start[0u]; nx < stop[0u]; ++nx) {
          kernel(Utils::Vector3i{nx, ny, nz}, linear_loop_index);
          linear_loop_index++;
        }
      }
    }
  }
}
