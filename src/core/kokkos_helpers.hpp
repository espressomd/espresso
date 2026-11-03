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

#include <algorithm>
#include <type_traits>

#include <Kokkos_Core.hpp>

/**
 * @brief Wrapper for ``Kokkos::deep_copy`` that skips fork/join when
 * the number of threads is 1.
 */
ESPRESSO_ATTR_ALWAYS_INLINE inline void
kokkos_deep_copy(auto const &exec_space, auto const &view, auto const &value) {
  using View = std::remove_cvref_t<decltype(view)>;
  if constexpr (Kokkos::SpaceAccessibility<
                    Kokkos::DefaultHostExecutionSpace,
                    typename View::memory_space>::accessible) {
    if (Kokkos::num_threads() == 1 and view.span_is_contiguous()) {
      std::fill_n(view.data(), view.span(),
                  static_cast<typename View::value_type>(value));
      return;
    }
  }
  Kokkos::deep_copy(exec_space, view, value);
}

/**
 * @brief Wrapper for ``Kokkos::parallel_for`` that skips fork/join when
 * the number of threads is 1.
 */
template <class ExecutionSpace = Kokkos::DefaultHostExecutionSpace>
ESPRESSO_ATTR_ALWAYS_INLINE inline void
kokkos_parallel_range_for(auto const &name, auto start, auto end,
                          auto const &kernel) {
  if (Kokkos::num_threads() > 1) {
    Kokkos::RangePolicy<ExecutionSpace> policy(start, end);
    Kokkos::parallel_for(name, policy, kernel);
  } else {
    for (auto p_index = start; p_index < end; ++p_index) {
      kernel(p_index);
    }
  }
}
