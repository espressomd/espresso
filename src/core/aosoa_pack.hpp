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

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM

#include "cell_system/CellStructure.hpp"

#include <Kokkos_Core.hpp>

struct CellStructure::AoSoA_pack {
  using PositionViewType =
      Kokkos::View<double *[3], Kokkos::LayoutRight, Kokkos::HostSpace>;
  using ChargeViewType = Kokkos::View<double *, Kokkos::HostSpace>;
  using IdViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using TypeViewType = Kokkos::View<int *, Kokkos::HostSpace>;
  using IdToIndexViewType = Kokkos::View<int *, Kokkos::HostSpace>;

  PositionViewType position;
  ChargeViewType charge;
  IdViewType id;
  TypeViewType type;
  IdToIndexViewType id_to_index;

  AoSoA_pack() = default;

  AoSoA_pack(std::size_t num_particles) { resize(num_particles); }

  void resize(std::size_t num_particles, int max_id = -1) {
    if (position.extent(0) == 0) {
      // First allocation
      position = PositionViewType("position", num_particles);
      charge = ChargeViewType("charge", num_particles);
      id = IdViewType("id", num_particles);
      type = TypeViewType("type", num_particles);
      id_to_index = IdToIndexViewType("id_to_index", num_particles);
    } else {
      // Reallocation
      Kokkos::realloc(position, num_particles);
      Kokkos::realloc(charge, num_particles);
      Kokkos::realloc(id, num_particles);
      Kokkos::realloc(type, num_particles);
    }
  }
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
