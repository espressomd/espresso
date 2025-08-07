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

#include "config/config.hpp"

#ifdef SHARED_MEMORY_PARALLELISM

#include "cell_system/CellStructure.hpp"

#include "aosoa_pack.hpp"
#include "custom_verlet_list.hpp"
#include "forces_cabana.hpp"

#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>

#include <iterator>
#include <utility>

inline void write_particle(Particle const &p, int const &id,
                           AoSoA_pack &aosoa) {
  aosoa.id(id) = p.id();
  aosoa.charge(id) = p.q();
  aosoa.type(id) = p.type();
  auto const &pos = p.pos();
  for (int d = 0; d < 3; ++d) {
    aosoa.position(id, d) = pos[d];
  }
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void construct_verlet_list(
    CellStructure &cell_structure, auto const &verlet_criterion,
    Kokkos::View<int *> const &id_to_index, const int max_id) {
  auto const &cells =
      std::as_const(cell_structure).decomposition().local_cells();
  auto const distance_function = detail::MinimalImageDistance{
      std::as_const(cell_structure).decomposition().box()};
  auto &verlet_list = cell_structure.get_cabana_verlet_list();

  // implementation detail: max_id refers to the max local particle id,
  // but ghost particles from other ranks may have larger particle ids;
  // in addition, -1 is used as a sentinel value for particle ids

  auto intra_kernel = [&cells, &distance_function, &verlet_criterion,
                       &id_to_index, &verlet_list, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto it = local_particles.begin(); it != local_particles.end(); ++it) {
      auto const &p1 = *it;
      if (p1.id() <= max_id) {
        auto const ii = id_to_index(p1.id());
        if (ii >= 0) {
          // pairs in this cell
          for (auto jt = std::next(it); jt != local_particles.end(); ++jt) {
            if ((*jt).id() <= max_id) {
              if (verlet_criterion(p1, *jt, distance_function(p1, *jt))) {
                auto const jj = id_to_index((*jt).id());
                if (jj >= 0) {
                  verlet_list.addNeighborLB(ii, jj);
                }
              }
            }
          }
        }
      }
    }
  };

  auto inter_kernel = [&cells, &distance_function, &verlet_criterion,
                       &id_to_index, &verlet_list, max_id](const int i) {
    auto &local_particles = cells[i]->particles();
    for (auto it = local_particles.begin(); it != local_particles.end(); ++it) {
      auto const &p1 = *it;
      if (p1.id() <= max_id) {
        auto const ii = id_to_index(p1.id());
        if (ii >= 0) {
          // pairs with neighboring cells
          for (auto &neighbor : cells[i]->neighbors().red()) {
            for (auto const &p2 : neighbor->particles()) {
              if (p2.id() <= max_id) {
                if (verlet_criterion(p1, p2, distance_function(p1, p2))) {
                  auto const jj = id_to_index(p2.id());
                  if (jj >= 0) {
                    verlet_list.addNeighbor(ii, jj);
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  Kokkos::parallel_for("inter", cells.size(), intra_kernel);
  Kokkos::fence();

  Kokkos::parallel_for("intra", cells.size(), inter_kernel);
  Kokkos::fence();
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_cabana_state(CellStructure &cell_structure, auto const &verlet_criterion,
                    double const pair_cutoff) {

  int num_threads = execution_space().concurrency();
  auto const rebuild = cell_structure.get_rebuild_cabana_verlet_list() or
                       (not cell_structure.use_verlet_list);

  if (rebuild) {
    // If we have to rebuild, we need to count the particles
    cell_structure.set_index_map(); // parallelized index_map
    // Create essential variables for MD
    cell_structure.rebuild_local_properties(num_threads, pair_cutoff);
  } else {
    // If we do not rebuild we can use the saved map
    cell_structure.reset_local_properties();
  }
  auto const unique_particles = cell_structure.get_unique_particles();
  auto aosoa = cell_structure.get_aosoa_data();
  auto max_id = cell_structure.get_cached_max_local_particle_id();

  // Fill the essential variables for MD
  {
    // ===================================================
    // Fill particle storage
    // ===================================================
    Kokkos::View<int *> id_to_index(
        Kokkos::ViewAllocateWithoutInitializing("id_to_index"), max_id + 1);
    Kokkos::deep_copy(id_to_index, -1);

    using policy_type = Kokkos::RangePolicy<execution_space>;
    Kokkos::parallel_for(
        "AoSoA write", policy_type(0, unique_particles.size()),
        [&unique_particles, &aosoa, &id_to_index](const int p_id) {
          write_particle(*unique_particles.at(p_id), p_id, aosoa);
          id_to_index(unique_particles.at(p_id)->id()) = p_id;
        });
    Kokkos::fence();

    // ===================================================
    // Get Verlet Pairs and Fill Verlet list
    // ===================================================

    if (rebuild) {
      construct_verlet_list(cell_structure, verlet_criterion, id_to_index,
                            max_id);
      cell_structure.mark_rebuild_cabana_verlet_list_as_UpToDate();
    }
  }
}

void cabana_short_range(auto const &bond_kernel, auto const &forces_kernel,
                        CellStructure &cell_structure, double pair_cutoff,
                        double bond_cutoff) {
  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
    cell_structure.bond_loop(bond_kernel);
  }

  // Cabana short range loop
  if (pair_cutoff > 0.) {
    auto &cabana_verlet_list = cell_structure.get_cabana_verlet_list();
    Kokkos::RangePolicy<execution_space> policy(
        0, cell_structure.get_unique_particles().size());
    Cabana::neighbor_parallel_for(policy, forces_kernel, cabana_verlet_list,
                                  Cabana::FirstNeighborsTag(),
                                  Cabana::SerialOpTag());
    Kokkos::fence();
  }
}

#endif // SHARED_MEMORY_PARALLELISM
