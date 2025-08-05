/*
 * Copyright (C) 2010-2024 The ESPResSo project
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

#include "Cell.hpp"
#include "config/config.hpp"

#include <span>
#include <vector>

#ifdef SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#endif

// Forward declaration
struct CellStructure;

/**
 * @brief Run a kernel on all local particles with enumeration.
 * The kernel is called with (index, particle) and is assumed to be thread-safe.
 *
 * @tparam Kernel Callable with signature void(int, Particle&)
 * @param cs The CellStructure containing the particles
 * @param kernel The kernel to apply to each particle with its index
 */
template <typename Kernel>
void enumerate_local_particles(CellStructure const &cs, Kernel &&kernel);

// Include the implementation
#include "CellStructure.hpp"

template <typename Kernel>
void enumerate_local_particles(CellStructure const &cs, Kernel &&kernel) {
#ifdef SHARED_MEMORY_PARALLELISM
  if (cs.use_parallel_for_each_local_particle()) {
    auto const local_cells = cs.decomposition().local_cells();

    // Step 1: Calculate cell offsets
    std::vector<int> cell_offsets(local_cells.size() + 1, 0);

    // Calculate cumulative sum of particles per cell
    for (size_t i = 0; i < local_cells.size(); ++i) {
      cell_offsets[i + 1] =
          cell_offsets[i] + local_cells[i]->particles().size();
    }

    // Step 2: Parallel loop over cells
    Kokkos::parallel_for(
        "enumerate_local_particles", local_cells.size(), [&](auto cell_idx) {
          auto const base_offset = cell_offsets[cell_idx];
          auto &cell_particles = local_cells[cell_idx]->particles();

          // Loop over particles in this cell
          for (size_t part_idx = 0; part_idx < cell_particles.size();
               ++part_idx) {
            int global_index = base_offset + part_idx;
            kernel(global_index, *(cell_particles.begin() + part_idx));
          }
        });
    return;
  }
#endif
  // Sequential fallback
  int index = 0;
  for (auto &p : cs.local_particles()) {
    kernel(index++, p);
  }
}