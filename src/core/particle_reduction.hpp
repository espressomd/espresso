/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

#include <concepts>
#include <utility>

namespace Reduction {

/** @brief Custom reduction in the form required by Kokkos */
template <typename ResultType, class Kernel, class Reduction>
class KokkosReducer {
public:
  // Kokkos reduction functors need the value_type typedef.
  // This is the type of the result of the reduction.
  using value_type = ResultType;

  // Just like with parallel_for functors, you may specify
  // an execution_space typedef. If not provided, Kokkos
  // will use the default execution space by default.

  // kernels to wrap
  Reduction reduction_op;
  Kernel kernel;
  KokkosReducer(Kernel kernel, Reduction reduction_op)
      : reduction_op(std::move(reduction_op)), kernel(std::move(kernel)) {}

  inline void operator()(std::integral auto const i, value_type &update) const {
    kernel(i, update);
  }

  // "Join" intermediate results from different threads.
  // This should normally implement the same reduction
  // operation as operator() above.
  inline void join(value_type &dst, value_type const &src) const {
    reduction_op(dst, src);
  }
};

template <typename ResultType, class Kernel, class Reduction>
KokkosReducer<ResultType, Kernel, Reduction>
make_kokkos_reducer(Kernel k, Reduction reduce_op) {
  return KokkosReducer<ResultType, Kernel, Reduction>(k, reduce_op);
}

} // namespace Reduction

/** @brief Run a reduction over all particles.
 *
 * @param cs  cell structure to iterate over
 * @param add_partial function that accumulates a result from a single particle
 * @param reduce_op function that joins two reduction results
 *
 * both functions have to implement the same reduction.
 */
template <typename ResultType>
ResultType reduce_over_local_particles(
    CellStructure const &cs,
    std::invocable<ResultType &, Particle const &> auto add_partial,
    [[maybe_unused]] std::invocable<ResultType &, ResultType const &> auto
        reduce_op) {

  static_assert(std::is_invocable_r_v<void, decltype(add_partial), ResultType &,
                                      Particle const &>);
  static_assert(std::is_invocable_r_v<void, decltype(reduce_op), ResultType &,
                                      ResultType const &>);
  ResultType result{};

  auto const &cells = cs.decomposition().local_cells();
  if (cells.size() > 1) { // parallel loop over cells
    auto reducer = Reduction::make_kokkos_reducer<ResultType>(
        [&cells, add_partial](std::size_t const c_index, ResultType &res) {
          for (auto const &p : cells[c_index]->particles()) {
            add_partial(res, p);
          }
        },
        reduce_op);
    Kokkos::parallel_reduce( // loop over cells
        "reduce_on_local_particle",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0},
                                                               cells.size()),
        reducer, result);
    return result;
  }
  // single cell case
  auto const &particles = cells.front()->particles();
  auto reducer = Reduction::make_kokkos_reducer<ResultType>(
      [&particles, add_partial](std::size_t const p_index, ResultType &res) {
        add_partial(res, std::as_const(*(particles.begin() + p_index)));
      },
      reduce_op);
  Kokkos::parallel_reduce( // loop over particles
      "reduce_on_local_particle",
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0},
                                                             particles.size()),
      reducer, result);
  return result;
}
