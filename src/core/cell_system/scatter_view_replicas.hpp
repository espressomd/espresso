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

/**
 * @file
 * Fused reduction and zeroing of the per-thread replicas of a duplicated
 * Kokkos ScatterView. Relies on Kokkos internals (the protected replica
 * member and its layout); @c scatter_view_replicas_test checks that the
 * fused path is taken for the force scatter view of @ref CellStructure.
 */

#include "kokkos_helpers.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <algorithm>
#include <cstddef>
#include <type_traits>

namespace ScatterViewReplicaReduction {
/** @brief Duplicated (one replica per thread) sum ScatterView. */
template <class> struct is_duplicated_sum_scatter_view : std::false_type {};
template <class DataType, class Layout, class Device, class Contribution>
struct is_duplicated_sum_scatter_view<Kokkos::Experimental::ScatterView<
    DataType, Layout, Device, Kokkos::Experimental::ScatterSum,
    Kokkos::Experimental::ScatterDuplicated, Contribution>> : std::true_type {};

/**
 * @brief Per-thread replicas of a duplicated ScatterView, a protected
 * member for which Kokkos (5.0) has no public accessor.
 */
template <class ScatterView> struct ScatterViewReplicas : ScatterView {
  static auto const &get(ScatterView const &scatter_view) {
    return scatter_view.*(&ScatterViewReplicas::internal_view);
  }
};

/**
 * @brief Number of values per block of the replica reduction. Each block is
 * reduced replica by replica, so that only two contiguous streams (the
 * destination block and one replica block) are accessed at a time. Iterating
 * over all replicas for every value instead accesses n_replicas streams that
 * are n_values apart; when n_values * sizeof(double) is a multiple of the
 * page size (for forces: a multiple of 512 particles), all streams map to
 * the same L1 cache sets and their loads and stores 4K-alias, which made
 * such a per-value loop 2.5x slower than contribute() + reset() at 32
 * threads and 320000 particles.
 */
inline constexpr std::size_t replica_reduction_block = 512u;

/**
 * @brief Add the replicas of a scatter view to its destination view and
 * zero them, in the summation order of the Kokkos reduction over the
 * replicas (dest += rep_0; dest += rep_1; ...).
 * Falls back to @c Kokkos::Experimental::contribute when the scatter view
 * has no replicas (non-duplicated, e.g. serial execution space), is not a
 * sum, or when the replicas are not stored as one contiguous array of
 * destination-sized blocks. A duplicated sum scatter view whose replicas or
 * destination are not LayoutRight is a compile-time error, not a fallback.
 * @return true if the replicas were zeroed.
 */
template <class ScatterView, class View>
bool contribute_and_reset_replicas(ScatterView const &scatter_view,
                                   View &dest) {
  if constexpr (is_duplicated_sum_scatter_view<ScatterView>::value) {
    using execution_space = typename ScatterView::execution_space;
    auto const &replicas = ScatterViewReplicas<ScatterView>::get(scatter_view);
    using ReplicasType = std::remove_cvref_t<decltype(replicas)>;
    static_assert(std::is_same_v<typename ReplicasType::array_layout,
                                 Kokkos::LayoutRight>);
    static_assert(
        std::is_same_v<typename View::array_layout, Kokkos::LayoutRight>);
    static_assert(ReplicasType::rank() == View::rank() + 1u);
    static_assert(std::is_same_v<typename ReplicasType::non_const_value_type,
                                 typename View::non_const_value_type>);
    auto const n_replicas = replicas.extent(0);
    auto const n_values = dest.size();
    auto layout_matches = replicas.span_is_contiguous() and
                          dest.span_is_contiguous() and
                          replicas.data() != dest.data() and
                          replicas.size() == n_replicas * n_values and
                          replicas.stride(0) == n_values;
    for (std::size_t r = 0u; r < View::rank(); ++r) {
      layout_matches &= (replicas.extent(r + 1u) == dest.extent(r));
    }
    if (layout_matches) {
      using value_type = typename View::non_const_value_type;
      auto *const rep = replicas.data();
      auto *const out = dest.data();
      constexpr auto block = replica_reduction_block;
      auto const n_blocks = (n_values + block - 1u) / block;
      kokkos_parallel_range_for<execution_space>(
          "ContributeResetDuplicates", std::size_t{0}, n_blocks,
          [rep, out, n_replicas, n_values](std::size_t const b) {
            auto const begin = b * block;
            auto const end = std::min(begin + block, n_values);
            for (std::size_t t = 0u; t < n_replicas; ++t) {
              auto *const rep_t = rep + t * n_values;
              // Both GCC and Clang auto-vectorize this loop (with a runtime
              // check that out and rep_t do not overlap). A forced
              // vectorization hint such as omp simd is an error under Clang
              // with -ftrapping-math (ESPRESSO_BUILD_WITH_FPE), where strict
              // floating-point semantics disable the loop vectorizer.
              for (std::size_t i = begin; i < end; ++i) {
                out[i] += rep_t[i];
                rep_t[i] = value_type{0};
              }
            }
          });
      execution_space{}.fence();
      return true;
    }
  }
  Kokkos::Experimental::contribute(dest, scatter_view);
  return false;
}

/** @brief Check that the replicas of a scatter view are all zero. */
template <class ScatterView>
bool replicas_are_zero(ScatterView const &scatter_view) {
  if constexpr (is_duplicated_sum_scatter_view<ScatterView>::value) {
    auto const &replicas = ScatterViewReplicas<ScatterView>::get(scatter_view);
    auto const *const data = replicas.data();
    return std::all_of(data, data + replicas.size(),
                       [](auto const value) { return value == 0; });
  }
  return true;
}
} // namespace ScatterViewReplicaReduction
