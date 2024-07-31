/*
 * Copyright (C) 2024 The ESPResSo project
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

#include <concepts>
#include <cstddef>

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
