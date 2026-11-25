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

#pragma once

#include <config/config.hpp>

#include <Kokkos_Core.hpp>

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

/**
 * @brief Flat list of the bonds between @p N particles that are evaluated
 * on this MPI rank.
 *
 * Each row holds the particles of one bond, the owner first. The rows hold
 * particle ids while the list is filled with @ref append(), and indices
 * into the particle arrays once @ref convert_ids_to_indices() was called.
 */
template <std::size_t N> struct LocalBondList {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  using ListType = Kokkos::View<int *[N], Kokkos::LayoutRight, execution_space>;
  using IdsType = Kokkos::View<int *, Kokkos::LayoutRight, execution_space>;
  using IdToIndex = Kokkos::View<int *, execution_space>;

  /** @brief Number of bonds, i.e. of rows that are filled. */
  int count = 0;
  /** @brief Per bond, its particles. */
  ListType list;
  /** @brief Per bond, the key of the bonded interaction. */
  IdsType ids;

  /** @brief Empty the list and make room for @p capacity bonds. */
  void allocate(int capacity);

  /** @brief Empty the list and release its memory. */
  void clear();

  /** @brief Empty the list, keeping its memory. */
  void reset();

  /**
   * @brief Append a bond, given as particle ids. Safe to call concurrently.
   * The capacity of the list must not be exceeded.
   */
  void append(int bond_id, std::span<const int> particle_ids) {
    assert(particle_ids.size() == N);
    auto const row = Kokkos::atomic_fetch_add(&count, 1);
    assert(static_cast<std::size_t>(row) < list.extent(0));
    for (std::size_t col = 0u; col < N; ++col) {
      list(row, col) = particle_ids[col];
    }
    ids(row) = bond_id;
  }

  /** @brief Replace the particle ids in all rows by particle indices. */
  void convert_ids_to_indices(IdToIndex const &id_to_index);

  /**
   * @brief Order the rows by their first entry, keeping the order of rows
   * with the same first entry. Applied to particle indices, the bond loops
   * then visit the particles in memory order.
   */
  void sort_by_first_particle();

#ifdef ESPRESSO_COLLISION_DETECTION
  /**
   * @brief Register a bond created after the list was filled.
   * It is counted right away, but only stored by @ref merge_new_bonds().
   */
  void add_new_bond(int bond_id, std::span<const int> particle_ids,
                    IdToIndex const &id_to_index);

  /** @brief Forget the bonds that were registered but not merged. */
  void discard_new_bonds();

  /** @brief Store the registered bonds at the end of the list. */
  void merge_new_bonds();

private:
  /** @brief Particle indices and bond ids of the registered bonds. */
  std::vector<int> m_new_list, m_new_ids;
#endif // ESPRESSO_COLLISION_DETECTION
};

/**
 * @brief Flat lists of the bonds evaluated on this MPI rank, by number of
 * particles.
 */
struct LocalBondState {
  using IdToIndex = LocalBondList<2>::IdToIndex;

  LocalBondList<2> pairs;
  LocalBondList<3> angles;
  LocalBondList<4> dihedrals;

  bool empty() const {
    return pairs.count == 0 and angles.count == 0 and dihedrals.count == 0;
  }

  /** @brief Empty the lists and make room for the given numbers of bonds. */
  void allocate(int n_pairs, int n_angles, int n_dihedrals);

  /** @brief Empty the lists and release their memory. */
  void clear();

  /** @brief Empty the lists, keeping their memory. */
  void reset();

  /**
   * @brief Append a bond, given as particle ids, to the list for its
   * number of particles. Safe to call concurrently. Bonds without
   * partners are not stored.
   */
  void append(int bond_id, std::span<const int> particle_ids) {
    switch (particle_ids.size()) {
    case 2u:
      pairs.append(bond_id, particle_ids);
      break;
    case 3u:
      angles.append(bond_id, particle_ids);
      break;
    case 4u:
      dihedrals.append(bond_id, particle_ids);
      break;
    default:
      break;
    }
  }

  /** @brief Replace the particle ids in all lists by particle indices. */
  void convert_ids_to_indices(IdToIndex const &id_to_index);

  /** @copydoc LocalBondList::sort_by_first_particle */
  void sort_by_first_particle();

#ifdef ESPRESSO_COLLISION_DETECTION
  /** @copydoc LocalBondList::add_new_bond */
  void add_new_bond(int bond_id, std::span<const int> particle_ids,
                    IdToIndex const &id_to_index);
  /** @copydoc LocalBondList::discard_new_bonds */
  void discard_new_bonds();
  /** @copydoc LocalBondList::merge_new_bonds */
  void merge_new_bonds();
#endif // ESPRESSO_COLLISION_DETECTION
};
