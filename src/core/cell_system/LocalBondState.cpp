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

#include <config/config.hpp>

#include "LocalBondState.hpp"

#include <Kokkos_Core.hpp>

#include <cassert>
#include <cstddef>
#include <span>
#include <utility>

template <std::size_t N> void LocalBondList<N>::allocate(int capacity) {
  // also drops pending new bonds, which would otherwise make the count
  // inconsistent in a later discard_new_bonds() or merge_new_bonds()
  reset();
  auto const n_rows = static_cast<std::size_t>(capacity);
  if (list.is_allocated()) {
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), list,
                    n_rows);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), ids,
                    n_rows);
  } else {
    list =
        ListType(Kokkos::view_alloc(execution_space{},
                                    Kokkos::WithoutInitializing, "bond_list"),
                 n_rows);
    ids = IdsType(Kokkos::view_alloc(execution_space{},
                                     Kokkos::WithoutInitializing, "bond_ids"),
                  n_rows);
  }
}

template <std::size_t N> void LocalBondList<N>::clear() {
  reset();
  // Reset Kokkos Views to default (unallocated) state
  list = ListType();
  ids = IdsType();
}

template <std::size_t N> void LocalBondList<N>::reset() {
  count = 0;
#ifdef ESPRESSO_COLLISION_DETECTION
  m_new_list.clear();
  m_new_ids.clear();
#endif
}

template <std::size_t N>
void LocalBondList<N>::convert_ids_to_indices(IdToIndex const &id_to_index) {
  if (count == 0) {
    return;
  }
  Kokkos::parallel_for("convert_bond_ids_to_indices",
                       Kokkos::RangePolicy<execution_space>(0, count),
                       [&bond_list = list, &id_to_index](int row) {
                         for (std::size_t col = 0u; col < N; ++col) {
                           bond_list(row, col) =
                               id_to_index(bond_list(row, col));
                         }
                       });
}

#ifdef ESPRESSO_COLLISION_DETECTION
template <std::size_t N>
void LocalBondList<N>::add_new_bond(int bond_id,
                                    std::span<const int> particle_ids,
                                    IdToIndex const &id_to_index) {
  assert(particle_ids.size() == N);
  for (auto const pid : particle_ids) {
    assert(pid >= 0 and static_cast<std::size_t>(pid) < id_to_index.extent(0));
    // The bounds are checked above. Reading the raw data keeps the
    // message formatting of the Kokkos bounds check out of the paths of
    // the static analyzer, which reports a false positive there.
    auto const index = id_to_index.data()[pid];
    assert(index != -1);
    m_new_list.emplace_back(index);
  }
  m_new_ids.emplace_back(bond_id);
  count++;
}

template <std::size_t N> void LocalBondList<N>::discard_new_bonds() {
  count -= static_cast<int>(m_new_ids.size());
  m_new_list.clear();
  m_new_ids.clear();
}

template <std::size_t N> void LocalBondList<N>::merge_new_bonds() {
  if (m_new_ids.empty()) {
    return;
  }
  auto const n_new = m_new_ids.size();
  auto const old_count = count - static_cast<int>(n_new);

  auto const n_rows = static_cast<std::size_t>(count);
  ListType merged_list(Kokkos::view_alloc(execution_space{},
                                          Kokkos::WithoutInitializing,
                                          "bond_list"),
                       n_rows);
  IdsType merged_ids(Kokkos::view_alloc(execution_space{},
                                        Kokkos::WithoutInitializing,
                                        "bond_ids"),
                     n_rows);

  Kokkos::deep_copy(
      Kokkos::subview(merged_list, std::make_pair(0, old_count), Kokkos::ALL()),
      Kokkos::subview(list, std::make_pair(0, old_count), Kokkos::ALL()));
  Kokkos::deep_copy(Kokkos::subview(merged_ids, std::make_pair(0, old_count)),
                    Kokkos::subview(ids, std::make_pair(0, old_count)));

  Kokkos::parallel_for(
      "merge_new_bonds",
      Kokkos::RangePolicy<execution_space>(std::size_t{0}, n_new),
      [&merged_list, &merged_ids, old_count, &new_list = m_new_list,
       &new_ids = m_new_ids](std::size_t i) {
        auto const row = old_count + static_cast<int>(i);
        for (std::size_t col = 0u; col < N; ++col) {
          merged_list(row, col) = new_list[i * N + col];
        }
        merged_ids(row) = new_ids[i];
      });

  list = merged_list;
  ids = merged_ids;
  m_new_list.clear();
  m_new_ids.clear();
}
#endif // ESPRESSO_COLLISION_DETECTION

template struct LocalBondList<2>;
template struct LocalBondList<3>;
template struct LocalBondList<4>;

void LocalBondState::allocate(int n_pairs, int n_angles, int n_dihedrals) {
  pairs.allocate(n_pairs);
  angles.allocate(n_angles);
  dihedrals.allocate(n_dihedrals);
}

void LocalBondState::clear() {
  pairs.clear();
  angles.clear();
  dihedrals.clear();
}

void LocalBondState::reset() {
  pairs.reset();
  angles.reset();
  dihedrals.reset();
}

void LocalBondState::convert_ids_to_indices(IdToIndex const &id_to_index) {
  if (empty()) {
    return;
  }
  pairs.convert_ids_to_indices(id_to_index);
  angles.convert_ids_to_indices(id_to_index);
  dihedrals.convert_ids_to_indices(id_to_index);
  Kokkos::fence();
}

#ifdef ESPRESSO_COLLISION_DETECTION
void LocalBondState::add_new_bond(int bond_id,
                                  std::span<const int> particle_ids,
                                  IdToIndex const &id_to_index) {
  switch (particle_ids.size()) {
  case 2u:
    pairs.add_new_bond(bond_id, particle_ids, id_to_index);
    break;
  case 3u:
    angles.add_new_bond(bond_id, particle_ids, id_to_index);
    break;
  case 4u:
    dihedrals.add_new_bond(bond_id, particle_ids, id_to_index);
    break;
  default:
    break;
  }
}

void LocalBondState::discard_new_bonds() {
  pairs.discard_new_bonds();
  angles.discard_new_bonds();
  dihedrals.discard_new_bonds();
}

void LocalBondState::merge_new_bonds() {
  pairs.merge_new_bonds();
  angles.merge_new_bonds();
  dihedrals.merge_new_bonds();
}
#endif // ESPRESSO_COLLISION_DETECTION
