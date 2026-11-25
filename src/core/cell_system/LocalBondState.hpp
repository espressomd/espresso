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
 * Each row holds the particles of one bond, the owner first, as indices
 * into the particle arrays.
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

  /**
   * @brief Replace the content of the list.
   * @param rows      Particle indices, @p N per bond, one bond after the
   *                  other.
   * @param bond_ids  Per bond, the key of the bonded interaction.
   */
  void assign(std::span<const int> rows, std::span<const int> bond_ids);

  /** @brief Empty the list and release its memory. */
  void clear();

  /** @brief Empty the list, keeping its memory. */
  void reset();

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
#endif // ESPRESSO_COLLISION_DETECTION

private:
  /** @brief Empty the list and make room for @p capacity bonds. */
  void allocate(int capacity);

#ifdef ESPRESSO_COLLISION_DETECTION
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

  /** @brief Empty the lists and release their memory. */
  void clear();

  /** @brief Empty the lists, keeping their memory. */
  void reset();

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
