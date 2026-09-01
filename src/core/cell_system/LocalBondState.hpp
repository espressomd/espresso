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

  // Per-particle gather structure for pair bonds: row = a local particle's
  // AoSoA index, one slot per BondList entry (primary or mirror, unlike the
  // per-bond lists below) it holds for a pair bond. Sized exactly to the
  // per-particle degree found during the counting pass (see
  // CellStructure::set_index_map), so no overflow handling is needed.
  // Columns: 0 = other participant id (until resolved to an AoSoA index),
  // 1 = bond_id, 2 = 1 if this is the primary (owning) entry, 0 if a mirror.
  using PPPairSlotType =
      Kokkos::View<int **[3], Kokkos::LayoutRight, execution_space>;
  using PPPairDegreeType =
      Kokkos::View<int *, Kokkos::LayoutRight, execution_space>;

  // Same idea as PPDihedralSlotType below, for angle bonds: not every
  // 3-body force function sharing calc_bonded_three_body_force's dispatch
  // has angle_generic_force's f_left(a,b) == f_right(b,a) symmetry (e.g.
  // IBMTriel does not), so a mirror-holding arm cannot just treat itself as
  // the "left" participant and get the right answer -- every row instead
  // stores the full original (vertex, arm1, arm2) plus which position
  // "self" occupies, with a mirror's position derived once at build time
  // (CellStructure::update_bond_storage) by looking up the owner's primary
  // entry, mirroring the dihedral chain_slot derivation. Columns: 0/1/2 =
  // (vertex, arm1, arm2) participant ids (until resolved to AoSoA
  // indices), 3 = bond_id, 4 = this row's own position (0=vertex, 1=arm1,
  // 2=arm2), 5 = bond_index: this bond's row in angle_list/angle_ids if its
  // primary entry is resolvable on this rank (i.e. the owning particle is
  // local, not a ghost, here -- see AngleBondsForceComputeKernel, which
  // already applies this row's share directly via a ScatterView in that
  // case), or -1 if it isn't (a bond straddling a rank boundary whose owner
  // lives elsewhere), in which case AngleBondsKernel falls back to
  // evaluating the bond directly from this row's own columns 0-4.
  using PPAngleSlotType =
      Kokkos::View<int **[6], Kokkos::LayoutRight, execution_space>;
  using PPAngleDegreeType =
      Kokkos::View<int *, Kokkos::LayoutRight, execution_space>;

  // Same idea for dihedral bonds. Unlike pair/angle, the dihedral force
  // formula gives each of the 4 chain positions an algebraically distinct
  // contribution (see calc_bonded_four_body_force), so there is no
  // self-referenced/symmetric shortcut: every row stores the full chain
  // (all 4 participant ids, in the bond's original creation order) plus
  // which position "self" occupies. For a mirror entry this position is
  // derived once at build time (CellStructure::update_bond_storage) by
  // looking up the owner's primary entry, since a mirror's own
  // partner_ids() alone cannot disambiguate which of the 3 non-owner
  // positions it originally was. Columns: 0/1/2/3 = the 4 chain participant
  // ids (until resolved to AoSoA indices), 4 = bond_id, 5 = this row's own
  // chain position (0..3), 6 = bond_index into dihedral_list/dihedral_ids,
  // or -1 -- same rank-boundary fallback rationale as PPAngleSlotType's
  // column 5.
  using PPDihedralSlotType =
      Kokkos::View<int **[7], Kokkos::LayoutRight, execution_space>;
  using PPDihedralDegreeType =
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

  // Number of rows (local particles) the per-particle structures below are
  // currently sized for.
  int pp_num_particles = 0;
  PPPairDegreeType pp_pair_degree;
  PPPairSlotType pp_pair_slots;
  PPAngleDegreeType pp_angle_degree;
  PPAngleSlotType pp_angle_slots;
  PPDihedralDegreeType pp_dihedral_degree;
  PPDihedralSlotType pp_dihedral_slots;

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

  /** Allocate or reallocate the per-particle pair-bond gather structure. */
  void allocate_pp_pair(int num_particles, int max_degree);

  /** Allocate or reallocate the per-particle angle-bond gather structure. */
  void allocate_pp_angle(int num_particles, int max_degree);

  /** Allocate or reallocate the per-particle dihedral-bond gather structure. */
  void allocate_pp_dihedral(int num_particles, int max_degree);

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
