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

#include "aosoa_pack.hpp"
#include "custom_verlet_list.hpp"
#include "forces_cabana.hpp"

#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>

#include <iterator>
#include <span>
#include <utility>

template <class KokkosRangePolicy = Kokkos::RangePolicy<>>
ESPRESSO_ATTR_ALWAYS_INLINE inline void
kokkos_parallel_range_for(auto const &name, auto start, auto end,
                          auto const &kernel) {
  if (Kokkos::num_threads() > 1) {
    KokkosRangePolicy policy(start, end);
    Kokkos::parallel_for(name, policy, kernel);
  } else {
    for (auto p_index = start; p_index < end; ++p_index) {
      kernel(p_index);
    }
  }
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void
commit_particle(Particle const &p, auto const index,
                CellStructure::AoSoA_pack &aosoa) {
  aosoa.id(index) = p.id();
#ifdef ESPRESSO_ELECTROSTATICS
  aosoa.charge(index) = p.q();
#endif
  aosoa.type(index) = p.type();
  auto const &pos = p.pos();
  aosoa.position(index, 0) = pos[0];
  aosoa.position(index, 1) = pos[1];
  aosoa.position(index, 2) = pos[2];
}

ESPRESSO_ATTR_ALWAYS_INLINE inline void construct_verlet_list(
    std::span<Cell *const> cells, BoxGeometry const &box_geo,
    CellStructure::ListType &verlet_list, auto const &verlet_criterion,
    Kokkos::View<int *> const &id_to_index, int const max_id) {

  auto const distance_function = detail::MinimalImageDistance{box_geo};

  // implementation detail: max_id refers to the max local particle id,
  // but ghost particles from other ranks may have larger particle ids;
  // -1 is used as a sentinel value for particle ids from other threads

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
  using execution_space = Kokkos::DefaultExecutionSpace;
  using policy_type = Kokkos::RangePolicy<execution_space>;
  auto const rebuild = cell_structure.prepare_verlet_list_cabana(pair_cutoff);
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto const max_id = cell_structure.get_cached_max_local_particle_id();
  auto &aosoa = cell_structure.get_aosoa();

  // ===================================================
  // Fill particle storage
  // ===================================================
  Kokkos::View<int *> id_to_index(
      Kokkos::ViewAllocateWithoutInitializing("id_to_index"), max_id + 1);
  Kokkos::deep_copy(id_to_index, -1);

  kokkos_parallel_range_for<policy_type>(
      "AoSoA write", std::size_t{0}, n_part,
      [&unique_particles, &aosoa, &id_to_index](int const index) {
        auto const &p = *unique_particles.at(index);
        commit_particle(p, index, aosoa);
        id_to_index(p.id()) = index;
      });
  Kokkos::fence();

  // ===================================================
  // Get Verlet pairs and fill Verlet list
  // ===================================================
  if (rebuild) {
    cell_structure.rebuild_verlet_list_cabana(
        [&](std::span<Cell *const> cells, BoxGeometry const &box,
            CellStructure::ListType &verlet_list) {
          construct_verlet_list(std::move(cells), box, verlet_list,
                                verlet_criterion, id_to_index, max_id);
        });
  }
}

#ifdef ESPRESSO_ELECTROSTATICS
ESPRESSO_ATTR_ALWAYS_INLINE inline void
update_aosoa_charges(CellStructure &cell_structure) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using policy_type = Kokkos::RangePolicy<execution_space>;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const n_part = unique_particles.size();
  auto &aosoa = cell_structure.get_aosoa();

  kokkos_parallel_range_for<policy_type>(
      "AoSoA update charges", std::size_t{0}, n_part,
      [&unique_particles, &aosoa](std::size_t const index) {
        aosoa.charge(index) = unique_particles.at(index)->q();
      });
}
#endif

void cabana_short_range(auto const &bond_kernel, auto const &forces_kernel,
                        CellStructure &cell_structure, double pair_cutoff,
                        double bond_cutoff) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
    cell_structure.bond_loop(bond_kernel);
  }

  // Cabana short range loop
  if (pair_cutoff > 0.) {
    auto const &verlet_list = cell_structure.get_verlet_list_cabana();
    Kokkos::RangePolicy<execution_space> policy(
        std::size_t{0}, cell_structure.get_unique_particles().size());
    Cabana::neighbor_parallel_for(policy, forces_kernel, verlet_list,
                                  Cabana::FirstNeighborsTag(),
                                  Cabana::SerialOpTag());
    Kokkos::fence();
  }
}

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
