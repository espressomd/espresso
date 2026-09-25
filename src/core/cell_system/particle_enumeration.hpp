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

#include "Cell.hpp"
#include "CellStructure.hpp"
#include "kokkos_helpers.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <numeric>
#include <vector>

/**
 * @brief Run a kernel on all local particles with enumeration.
 * The kernel is called with (index, particle) and is assumed to be thread-safe.
 *
 * @tparam Kernel Callable with signature <tt>void(std::size_t, Particle&)</tt>
 * @param cs The cell structure containing the particles
 * @param kernel The kernel to apply to each particle with its index
 */
template <typename Kernel>
inline void enumerate_local_particles(CellStructure const &cs,
                                      Kernel &&kernel) {
  if (cs.use_parallel_for_each_local_particle()) {
    auto const local_cells = cs.decomposition().local_cells();

    // Use count() (raw committed rows) as the basis for cell_offsets to match
    // the inner loop, which iterates [offset, offset+count(). The function's
    // contract requires a clean store (no pending-removed rows); on a clean
    // store count() == particles().size() (live count), so behavior is
    // identical at all call sites.
    std::vector<std::size_t> cell_offsets(local_cells.size(), std::size_t{0});
    std::exclusive_scan(
        local_cells.begin(), local_cells.end(), cell_offsets.begin(),
        std::size_t{0},
        [](auto acc, auto const &cell) { return acc + cell->count(); });

    kokkos_parallel_range_for(
        "enumerate_local_particles", std::size_t{0}, local_cells.size(),
        [&](auto cell_idx) {
          auto const base_offset = cell_offsets[cell_idx];
          // The cell's committed rows are the contiguous range
          // [offset, offset+count); index it directly (this runs on a clean
          // store, so no pending-removed rows). Each cell gets its own view so
          // this stays thread-safe under the parallel_for over cells.
          auto *cell = local_cells[cell_idx];
          auto &store = cell->store();
          auto const row_offset = cell->offset();
          auto const n_part = cell->count();
          // Reuse one cached view across this cell's rows (this inner loop is
          // sequential -- one thread per cell -- so a single rebound view is
          // safe), REBOUND per row via attach_to_store instead of constructing
          // a fresh Particle per row. The heap-owning members stay
          // default-constructed and are never read while attached.
          Particle view;
          for (std::size_t p_index{0}; p_index < n_part; ++p_index) {
            auto global_index = base_offset + p_index;
            view.attach_to_store(store, static_cast<int>(row_offset + p_index));
            kernel(global_index, view);
          }
        });
    return;
  }
  // Sequential fallback
  std::size_t index = 0;
  for (auto &p : cs.local_particles()) {
    kernel(index++, p);
  }
}
