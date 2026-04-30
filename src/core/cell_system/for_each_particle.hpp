
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

template <typename Callable>
inline void
CellStructure::parallel_for_each_particle_impl(std::span<Cell *const> cells,
                                               Callable &f) const {
  if (cells.size() > 1) {
    Kokkos::parallel_for( // loop over cells
        "for_each_local_particle", cells.size(), [&](auto cell_idx) {
          for (auto &p : cells[cell_idx]->particles())
            f(p);
        });
  } else if (cells.size() == 1) {
    auto &particles = cells.front()->particles();
    Kokkos::parallel_for( // loop over particles
        "for_each_local_particle", particles.size(),
        [&](auto part_idx) { f(*(particles.begin() + part_idx)); });
  }
}
