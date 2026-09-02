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

#include <vector>

struct LocalBondState {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  using PairBondlistType =
      Kokkos::View<int *[2], Kokkos::LayoutRight, execution_space>;
  using PairBondIDType =
      Kokkos::View<int *, Kokkos::LayoutRight, execution_space>;
  using AngleBondlistType =
      Kokkos::View<int *[3], Kokkos::LayoutRight, execution_space>;
  using AngleBondIDType =
      Kokkos::View<int *, Kokkos::LayoutRight, execution_space>;
  using DihedralBondlistType =
      Kokkos::View<int *[4], Kokkos::LayoutRight, execution_space>;
  using DihedralBondIDType =
      Kokkos::View<int *, Kokkos::LayoutRight, execution_space>;

  // Bond counts
  int pair_count = 0;
  int angle_count = 0;
  int dihedral_count = 0;

  // Kokkos bond lists
  PairBondlistType pair_list;
  PairBondIDType pair_ids;
  AngleBondlistType angle_list;
  AngleBondIDType angle_ids;
  DihedralBondlistType dihedral_list;
  DihedralBondIDType dihedral_ids;

#ifdef ESPRESSO_COLLISION_DETECTION
  std::vector<int> new_pair_list, new_pair_ids;
  std::vector<int> new_angle_list, new_angle_ids;
  std::vector<int> new_dihedral_list, new_dihedral_ids;
#endif

  void reset_counts() {
    pair_count = 0;
    angle_count = 0;
    dihedral_count = 0;
  }

  void set_counts(int p, int a, int d) {
    pair_count = p;
    angle_count = a;
    dihedral_count = d;
  }

  /** Allocate or reallocate all Kokkos Views to current counts. */
  void allocate();

  /** Deallocates Views */
  void clear();

  /** Reset counts + collision vectors */
  void reset();

#ifdef ESPRESSO_COLLISION_DETECTION
  void clear_new_bonds();
  void add_new_bond(int bond_id, std::vector<int> const &particle_ids,
                    Kokkos::View<int *, execution_space> const &id_to_index);
  void rebuild();
#endif
};
