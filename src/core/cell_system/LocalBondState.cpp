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
    pair_list = PairBondlistType(
        Kokkos::ViewAllocateWithoutInitializing("pair_bond_list"), pair_count);
    pair_ids = PairBondIDType(
        Kokkos::ViewAllocateWithoutInitializing("pair_bond_id"), pair_count);
    angle_list = AngleBondlistType(
        Kokkos::ViewAllocateWithoutInitializing("angle_bond_list"),
        angle_count);
    angle_ids = AngleBondIDType(
        Kokkos::ViewAllocateWithoutInitializing("angle_bond_id"), angle_count);
    dihedral_list = DihedralBondlistType(
        Kokkos::ViewAllocateWithoutInitializing("dihedral_bond_list"),
        dihedral_count);
    dihedral_ids = DihedralBondIDType(
        Kokkos::ViewAllocateWithoutInitializing("dihedral_bond_id"),
        dihedral_count);
  }
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

void LocalBondState::add_new_bond(int bond_id,
                                  std::vector<int> const &particle_ids,
                                  Kokkos::View<int *> const &id_to_index) {
  if (particle_ids.size() == 2u) {
    new_pair_list.reserve(new_pair_list.size() + 2u);
    for (auto pid : particle_ids)
      new_pair_list.emplace_back(id_to_index(pid));
    new_pair_ids.emplace_back(bond_id);
    pair_count++;
  } else if (particle_ids.size() == 3u) {
    new_angle_list.reserve(new_angle_list.size() + 3u);
    for (auto pid : particle_ids)
      new_angle_list.emplace_back(id_to_index(pid));
    new_angle_ids.emplace_back(bond_id);
    angle_count++;
  } else if (particle_ids.size() == 4u) {
    new_dihedral_list.reserve(new_dihedral_list.size() + 4u);
    for (auto pid : particle_ids)
      new_dihedral_list.emplace_back(id_to_index(pid));
    new_dihedral_ids.emplace_back(bond_id);
    dihedral_count++;
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

  BondListT rebuilt_list(
      Kokkos::ViewAllocateWithoutInitializing("bond_list_rebuild"),
      total_bond_count);
  BondIDT rebuilt_ids(
      Kokkos::ViewAllocateWithoutInitializing("bond_id_rebuild"),
      total_bond_count);

  Kokkos::deep_copy(
      Kokkos::subview(rebuilt_list, std::make_pair(0, old_count),
                      Kokkos::ALL()),
      Kokkos::subview(bond_list, std::make_pair(0, old_count), Kokkos::ALL()));
  Kokkos::deep_copy(Kokkos::subview(rebuilt_ids, std::make_pair(0, old_count)),
                    Kokkos::subview(bond_ids, std::make_pair(0, old_count)));

  Kokkos::parallel_for(
      "copy_bondlist", new_bond_list.size(),
      [&bond_view = rebuilt_list, old_count, &new_data_view](auto flat_idx) {
        constexpr int NCols =
            BondListT::rank == 2 ? static_cast<int>(BondListT::static_extent(1))
                                 : 1;
        auto bond_idx = old_count + static_cast<int>(flat_idx / NCols);
        auto col_idx = static_cast<int>(flat_idx % NCols);
        bond_view(bond_idx, col_idx) = new_data_view(flat_idx);
      });

  Kokkos::parallel_for(
      "copy_bond_ids", new_bond_ids.size(),
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
