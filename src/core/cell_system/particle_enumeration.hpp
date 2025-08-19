/*
 * Copyright (C) 2025 The ESPResSo project
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

#ifdef SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#endif

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
#ifdef SHARED_MEMORY_PARALLELISM
  if (cs.use_parallel_for_each_local_particle()) {
    auto const local_cells = cs.decomposition().local_cells();

    std::vector<std::size_t> cell_offsets(local_cells.size(), std::size_t{0});
    std::exclusive_scan(local_cells.begin(), local_cells.end(),
                        cell_offsets.begin(), std::size_t{0},
                        [](auto acc, auto const &cell) {
                          return acc + cell->particles().size();
                        });

    Kokkos::parallel_for(
        "enumerate_local_particles", local_cells.size(), [&](auto cell_idx) {
          auto const base_offset = cell_offsets[cell_idx];
          auto &cell_particles = local_cells[cell_idx]->particles();
          auto const n_part = cell_particles.size();
          for (std::size_t p_index{0}; p_index < n_part; ++p_index) {
            auto global_index = base_offset + p_index;
            kernel(global_index, *(cell_particles.begin() + p_index));
          }
        });
    return;
  }
#endif // SHARED_MEMORY_PARALLELISM
  // Sequential fallback
  std::size_t index = 0;
  for (auto &p : cs.local_particles()) {
    kernel(index++, p);
  }
}
