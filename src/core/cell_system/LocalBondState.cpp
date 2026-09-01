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

#include <array>

void LocalBondState::allocate() {
  if (pair_list.is_allocated()) {
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), pair_list,
                    pair_count);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), pair_ids,
                    pair_count);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), angle_list,
                    angle_count);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), angle_ids,
                    angle_count);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    dihedral_list, dihedral_count);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    dihedral_ids, dihedral_count);
  } else {
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    pair_list = PairBondlistType(Kokkos::view_alloc(execution_space{},
                                                    Kokkos::WithoutInitializing,
                                                    "pair_bond_list"),
                                 pair_count);
    pair_ids = PairBondIDType(Kokkos::view_alloc(execution_space{},
                                                 Kokkos::WithoutInitializing,
                                                 "pair_bond_id"),
                              pair_count);
    angle_list = AngleBondlistType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "angle_bond_list"),
        angle_count);
    angle_ids = AngleBondIDType(Kokkos::view_alloc(execution_space{},
                                                   Kokkos::WithoutInitializing,
                                                   "angle_bond_id"),
                                angle_count);
    dihedral_list = DihedralBondlistType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "dihedral_bond_list"),
        dihedral_count);
    dihedral_ids = DihedralBondIDType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "dihedral_bond_id"),
        dihedral_count);
  }
}

void LocalBondState::allocate_pp_pair(int num_particles, int max_degree) {
  if (pp_pair_slots.is_allocated()) {
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    pp_pair_degree, num_particles);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    pp_pair_slots, num_particles, max_degree);
  } else {
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    pp_pair_degree = PPPairDegreeType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "pp_pair_degree"),
        num_particles);
    pp_pair_slots = PPPairSlotType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "pp_pair_slots"),
        num_particles, max_degree);
  }
  pp_num_particles = num_particles;
}

void LocalBondState::allocate_pp_angle(int num_particles, int max_degree) {
  if (pp_angle_slots.is_allocated()) {
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    pp_angle_degree, num_particles);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    pp_angle_slots, num_particles, max_degree);
  } else {
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    pp_angle_degree = PPAngleDegreeType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "pp_angle_degree"),
        num_particles);
    pp_angle_slots = PPAngleSlotType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "pp_angle_slots"),
        num_particles, max_degree);
  }
  pp_num_particles = num_particles;
}

void LocalBondState::allocate_pp_dihedral(int num_particles, int max_degree) {
  if (pp_dihedral_slots.is_allocated()) {
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    pp_dihedral_degree, num_particles);
    Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    pp_dihedral_slots, num_particles, max_degree);
  } else {
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    pp_dihedral_degree = PPDihedralDegreeType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "pp_dihedral_degree"),
        num_particles);
    pp_dihedral_slots = PPDihedralSlotType(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "pp_dihedral_slots"),
        num_particles, max_degree);
  }
  pp_num_particles = num_particles;
}

void LocalBondState::clear() {
  reset_counts();
  // Reset Kokkos Views to default (unallocated) state
  pair_list = PairBondlistType();
  pair_ids = PairBondIDType();
  angle_list = AngleBondlistType();
  angle_ids = AngleBondIDType();
  dihedral_list = DihedralBondlistType();
  dihedral_ids = DihedralBondIDType();
  pp_pair_degree = PPPairDegreeType();
  pp_pair_slots = PPPairSlotType();
  pp_angle_degree = PPAngleDegreeType();
  pp_angle_slots = PPAngleSlotType();
  pp_dihedral_degree = PPDihedralDegreeType();
  pp_dihedral_slots = PPDihedralSlotType();
  pp_num_particles = 0;
#ifdef ESPRESSO_COLLISION_DETECTION
  clear_new_bonds();
#endif
}

void LocalBondState::reset() {
  reset_counts();
#ifdef ESPRESSO_COLLISION_DETECTION
  clear_new_bonds();
#endif
}

#ifdef ESPRESSO_COLLISION_DETECTION
void LocalBondState::clear_new_bonds() {
  new_pair_list.clear();
  new_pair_ids.clear();
  new_angle_list.clear();
  new_angle_ids.clear();
  new_dihedral_list.clear();
  new_dihedral_ids.clear();
}

namespace {
// Resolves @p pid to an AoSoA index, or -1 if @p pid is not known on this
// rank at all (neither local nor ghost) -- id_to_index is sized to the
// cached max local particle id and default-filled with -1, so an unknown
// id can also be out of its bounds entirely, not just map to -1. Also
// rejects (as -1) an index that id_to_index resolves validly but which is
// stale relative to the pp_* structures' own row count (@p num_rows,
// i.e. LocalBondState::pp_num_particles): id_to_index tracks whichever
// particles/ghosts are currently committed to the AoSoA, which can outgrow
// the pp_* structures between two full rebuilds if ghost communication
// runs in between -- in that window the hot-add fast path must decline
// rather than write out of bounds; the bond is still picked up correctly
// at the next full rebuild.
template <typename IdToIndexT>
int safe_resolve_index(IdToIndexT const &id_to_index, int pid, int num_rows) {
  if (pid < 0 or static_cast<std::size_t>(pid) >= id_to_index.extent(0)) {
    return -1;
  }
  auto const idx = id_to_index(pid);
  if (idx < 0 or idx >= num_rows) {
    return -1;
  }
  return idx;
}

// Appends one gather-row for participant @p self_idx to a pp_* structure,
// growing the degree (2nd) dimension for every particle via Kokkos::resize
// (which -- unlike Kokkos::realloc, used by the full-rebuild allocate_pp_*()
// -- preserves and correctly reindexes existing data) if this particle's
// row is already full. Mirrors, for the hot-add path, what a full rebuild's
// counting + allocate_pp_*() + update_bond_storage() pass would produce for
// this one new entry.
template <typename DegreeT, typename SlotT, typename... Cols>
void append_pp_row(DegreeT &degree_view, SlotT &slot_view, int self_idx,
                   Cols... cols) {
  auto const degree = degree_view(self_idx);
  if (degree >= static_cast<int>(slot_view.extent(1))) {
    Kokkos::resize(slot_view, slot_view.extent(0), degree + 1);
  }
  int col = 0;
  ((slot_view(self_idx, degree, col++) = cols), ...);
  degree_view(self_idx) = degree + 1;
}
} // anonymous namespace

void LocalBondState::add_new_bond(
    int bond_id, std::vector<int> const &particle_ids,
    Kokkos::View<int *, execution_space> const &id_to_index) {
  if (particle_ids.size() == 2u) {
    new_pair_list.reserve(new_pair_list.size() + 2u);
    for (auto pid : particle_ids)
      new_pair_list.emplace_back(
          safe_resolve_index(id_to_index, pid, pp_num_particles));
    new_pair_ids.emplace_back(bond_id);
    pair_count++;

    // Every pp_pair row for this bond needs BOTH participants resolved to
    // a valid local/ghost AoSoA index (the "other" column is used directly
    // as an index by the force kernel) -- if either is unknown on this
    // rank, skip the gather fast-path for this bond entirely here; it will
    // be picked up correctly at the next full rebuild instead, exactly as
    // for a bond whose partner isn't yet ghost-visible via ::add_bond.
    auto const idx0 =
        safe_resolve_index(id_to_index, particle_ids[0], pp_num_particles);
    auto const idx1 =
        safe_resolve_index(id_to_index, particle_ids[1], pp_num_particles);
    if (idx0 >= 0 and idx1 >= 0) {
      // Column 3 (bond_index): safe_resolve_index() above already rejected
      // either participant being a ghost (out of [0, pp_num_particles)), so
      // both rows are guaranteed real/local -- PairBondsForceComputeKernel
      // will therefore always cover both shares once rebuild_bond_list()
      // (called once per step, after all of this step's collisions/hot-adds
      // are processed) merges this bond into pair_list. pair_count was just
      // incremented above, so pair_count - 1 is exactly this bond's future
      // row there -- known now, without waiting for rebuild_bond_list() to
      // determine it. Leaving this column unset here (as an earlier version
      // did) left it holding whatever the last full rebuild's
      // Kokkos::realloc(..., WithoutInitializing, ...) happened to leave in
      // that not-yet-written slot -- interpreted by the gather kernels as
      // "already resolved, skip" only by chance, and as "skip" is the
      // deliberately correct value here regardless, but as "unresolved,
      // fall back and re-evaluate" whenever that leftover happened to be
      // negative, double-counting this bond's force for however many steps
      // remained until the next full rebuild silently overwrote it.
      auto const bond_list_index = pair_count - 1;
      append_pp_row(pp_pair_degree, pp_pair_slots, idx0, idx1, bond_id, 1,
                    bond_list_index);
      append_pp_row(pp_pair_degree, pp_pair_slots, idx1, idx0, bond_id, 0,
                    bond_list_index);
    }
  } else if (particle_ids.size() == 3u) {
    new_angle_list.reserve(new_angle_list.size() + 3u);
    for (auto pid : particle_ids)
      new_angle_list.emplace_back(
          safe_resolve_index(id_to_index, pid, pp_num_particles));
    new_angle_ids.emplace_back(bond_id);
    angle_count++;

    auto const vertex_idx =
        safe_resolve_index(id_to_index, particle_ids[0], pp_num_particles);
    auto const arm1_idx =
        safe_resolve_index(id_to_index, particle_ids[1], pp_num_particles);
    auto const arm2_idx =
        safe_resolve_index(id_to_index, particle_ids[2], pp_num_particles);
    if (vertex_idx >= 0 and arm1_idx >= 0 and arm2_idx >= 0) {
      // Column 5 (bond_index): see the pair-bond case above for the
      // rationale and the bug this avoids -- same idea, all three
      // participants guaranteed real/local by the checks above, so
      // AngleBondsForceComputeKernel always covers all three shares once
      // rebuild_bond_list() merges this bond into angle_list, at the row
      // index angle_count - 1 (already known here, angle_count having just
      // been incremented above).
      auto const bond_list_index = angle_count - 1;
      append_pp_row(pp_angle_degree, pp_angle_slots, vertex_idx, vertex_idx,
                    arm1_idx, arm2_idx, bond_id, 0, bond_list_index);
      append_pp_row(pp_angle_degree, pp_angle_slots, arm1_idx, vertex_idx,
                    arm1_idx, arm2_idx, bond_id, 1, bond_list_index);
      append_pp_row(pp_angle_degree, pp_angle_slots, arm2_idx, vertex_idx,
                    arm1_idx, arm2_idx, bond_id, 2, bond_list_index);
    }
  } else if (particle_ids.size() == 4u) {
    new_dihedral_list.reserve(new_dihedral_list.size() + 4u);
    for (auto pid : particle_ids)
      new_dihedral_list.emplace_back(
          safe_resolve_index(id_to_index, pid, pp_num_particles));
    new_dihedral_ids.emplace_back(bond_id);
    dihedral_count++;

    std::array<int, 4> chain_idx{};
    bool all_valid = true;
    for (int c = 0; c < 4; ++c) {
      chain_idx[c] =
          safe_resolve_index(id_to_index, particle_ids[c], pp_num_particles);
      all_valid = all_valid and chain_idx[c] >= 0;
    }
    if (all_valid) {
      // Column 6 (bond_index): see the pair-bond case above for the
      // rationale. All 4 chain positions guaranteed real/local by
      // all_valid, so DihedralBondsForceComputeKernel always covers every
      // share once rebuild_bond_list() merges this bond into
      // dihedral_list, at the row index dihedral_count - 1.
      auto const bond_list_index = dihedral_count - 1;
      for (int c = 0; c < 4; ++c) {
        append_pp_row(pp_dihedral_degree, pp_dihedral_slots, chain_idx[c],
                      chain_idx[0], chain_idx[1], chain_idx[2], chain_idx[3],
                      bond_id, c, bond_list_index);
      }
    }
  }
}

namespace {
template <typename BondListT, typename BondIDT>
void rebuild_bond_list_impl(std::vector<int> const &new_bond_list,
                            std::vector<int> const &new_bond_ids,
                            BondListT &bond_list, BondIDT &bond_ids,
                            int total_bond_count) {
  if (new_bond_list.empty())
    return;

  auto new_data_view = Kokkos::View<const int *, Kokkos::HostSpace,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      new_bond_list.data(), new_bond_list.size());
  auto new_id_view = Kokkos::View<const int *, Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      new_bond_ids.data(), new_bond_ids.size());

  auto old_count = total_bond_count - static_cast<int>(new_bond_ids.size());

  using execution_space = Kokkos::DefaultHostExecutionSpace;
  BondListT rebuilt_list(Kokkos::view_alloc(execution_space{},
                                            Kokkos::WithoutInitializing,
                                            "bond_list_rebuild"),
                         total_bond_count);
  BondIDT rebuilt_ids(Kokkos::view_alloc(execution_space{},
                                         Kokkos::WithoutInitializing,
                                         "bond_id_rebuild"),
                      total_bond_count);

  Kokkos::deep_copy(
      Kokkos::subview(rebuilt_list, std::make_pair(0, old_count),
                      Kokkos::ALL()),
      Kokkos::subview(bond_list, std::make_pair(0, old_count), Kokkos::ALL()));
  Kokkos::deep_copy(Kokkos::subview(rebuilt_ids, std::make_pair(0, old_count)),
                    Kokkos::subview(bond_ids, std::make_pair(0, old_count)));

  Kokkos::parallel_for(
      "copy_bondlist",
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
          std::size_t{0}, new_bond_list.size()),
      [&bond_view = rebuilt_list, old_count, &new_data_view](auto flat_idx) {
        constexpr int NCols =
            BondListT::rank == 2 ? static_cast<int>(BondListT::static_extent(1))
                                 : 1;
        auto bond_idx = old_count + static_cast<int>(flat_idx / NCols);
        auto col_idx = static_cast<int>(flat_idx % NCols);
        bond_view(bond_idx, col_idx) = new_data_view(flat_idx);
      });

  Kokkos::parallel_for(
      "copy_bond_ids",
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
          std::size_t{0}, new_bond_ids.size()),
      [&id_view = rebuilt_ids, old_count, &new_id_view](auto idx) {
        id_view(old_count + static_cast<int>(idx)) = new_id_view(idx);
      });

  bond_list = rebuilt_list;
  bond_ids = rebuilt_ids;
}
} // anonymous namespace

void LocalBondState::rebuild() {
  rebuild_bond_list_impl(new_pair_list, new_pair_ids, pair_list, pair_ids,
                         pair_count);
  rebuild_bond_list_impl(new_angle_list, new_angle_ids, angle_list, angle_ids,
                         angle_count);
  rebuild_bond_list_impl(new_dihedral_list, new_dihedral_ids, dihedral_list,
                         dihedral_ids, dihedral_count);
  clear_new_bonds();
}
#endif // ESPRESSO_COLLISION_DETECTION
