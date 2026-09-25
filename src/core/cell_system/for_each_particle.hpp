
/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
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

#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

#include <span>

inline void CellStructure::parallel_for_each_particle_impl(
    std::span<Cell *const> cells, ParticleCallback auto &&f) const {
  auto &store = const_cast<ParticleStore &>(m_particle_store);
  if (cells.size() > 1) {
    Kokkos::parallel_for( // loop over cells
        "for_each_local_particle",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0},
                                                               cells.size()),
        [&](auto cell_idx) {
          // One reused view per cell (per thread), REBOUND per row via
          // attach_to_store instead of relying on the row-range iterator to
          // materialise a Particle per cell. The cell's committed rows are the
          // contiguous range [offset, offset+count); this runs on a clean store
          // (no pending-removed rows), so index it directly. The heap-owning
          // members stay default-constructed and are never read while attached.
          auto const offset = cells[cell_idx]->offset();
          auto const n_part = cells[cell_idx]->count();
          Particle p;
          for (std::size_t idx = 0u; idx < n_part; ++idx) {
            p.attach_to_store(store, static_cast<int>(offset + idx));
            f(p);
          }
        });
  } else if (cells.size() == 1) {
    // Single-cell parallel-over-particles: each index gets its OWN view (a
    // shared cached-view iterator would not be thread-safe here).
    auto const offset = cells.front()->offset();
    Kokkos::parallel_for( // loop over particles
        "for_each_local_particle",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
            std::size_t{0}, cells.front()->count()),
        [&](auto part_idx) {
          auto view = store.make_view(static_cast<int>(offset) + part_idx);
          f(view);
        });
  }
}

template <typename RowKernel>
inline void
CellStructure::parallel_for_each_local_row_impl(std::span<Cell *const> cells,
                                                RowKernel &kernel) const {
  // Column-kernel launcher: same iteration structure as
  // parallel_for_each_particle_impl (above) but hands the kernel the raw STORE
  // ROW instead of a rebound Particle view. The kernel body reads its hoisted
  // *_view() column handles directly by row -- no per-element view_host() /
  // address / stride recompute. The row set and traversal order match the view
  // path exactly (local cells tile [0, n_local) contiguously in cell order).
  if (cells.size() > 1) {
    Kokkos::parallel_for( // loop over cells
        "for_each_local_particle_row",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(std::size_t{0},
                                                               cells.size()),
        [&](auto cell_idx) {
          auto const offset = cells[cell_idx]->offset();
          auto const n_part = cells[cell_idx]->count();
          for (std::size_t idx = 0u; idx < n_part; ++idx) {
            kernel(static_cast<int>(offset + idx));
          }
        });
  } else if (cells.size() == 1) {
    auto const offset = cells.front()->offset();
    Kokkos::parallel_for( // loop over particles
        "for_each_local_particle_row",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
            std::size_t{0}, cells.front()->count()),
        [&](auto part_idx) {
          kernel(static_cast<int>(offset) + static_cast<int>(part_idx));
        });
  }
}
