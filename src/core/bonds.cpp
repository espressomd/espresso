/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

/**
 * @brief Ids of all participants of a bond except the one at @p excluded,
 * in their original relative order.
 */
static std::vector<int> other_participant_ids(std::vector<int> const &ids,
                                              std::size_t excluded) {
  std::vector<int> others;
  others.reserve(ids.size() - 1);
  for (std::size_t j = 0; j < ids.size(); ++j) {
    if (j != excluded) {
      others.push_back(ids[j]);
    }
  }
  return others;
}

/**
 * @brief Find a bond list entry matching @p bond_id whose partner ids
 * (order-independent) equal @p sorted_others, optionally restricted to a
 * specific role.
 */
static BondList::const_iterator
find_matching_bond(BondList const &bonds, int bond_id,
                   std::vector<int> const &sorted_others,
                   std::optional<bool> want_primary) {
  for (auto it = bonds.begin(); it != bonds.end(); ++it) {
    if (it->bond_id() != bond_id or
        it->partner_ids().size() != sorted_others.size()) {
      continue;
    }
    if (want_primary and it->is_primary() != *want_primary) {
      continue;
    }
    std::vector<int> candidate(it->partner_ids().begin(),
                               it->partner_ids().end());
    std::ranges::sort(candidate);
    if (candidate == sorted_others) {
      return it;
    }
  }
  return bonds.end();
}

bool add_bond(System::System &system, int bond_id,
              std::vector<int> const &particle_ids) {
  auto &cell_structure = *system.cell_structure;
  bool added = false;

  if (Particle *p = cell_structure.get_local_particle(particle_ids[0])) {
    // Primary entry: used for force/energy calculation, holds the other
    // participants' ids in creation order.
    BondView bond(bond_id, {particle_ids.data() + 1, particle_ids.size() - 1});
    p->bonds().insert(bond);
    added = true;
  }

  // Mirror entries, so the bond can be found and removed from any
  // participant.
  for (std::size_t i = 1; i < particle_ids.size(); ++i) {
    Particle *p = cell_structure.get_local_particle(particle_ids[i]);
    if (not p) {
      continue;
    }
    auto const others = other_participant_ids(particle_ids, i);
    BondView bond(bond_id, {others.data(), others.size()}, false);
    p->bonds().insert(bond);
    added = true;
  }

  return added;
}

bool remove_bond(System::System &system, int bond_id,
                 std::vector<int> const &particle_ids, int skip_id) {
  auto &cell_structure = *system.cell_structure;

  auto try_remove = [&](int pid, std::vector<int> const &sorted_others,
                        std::optional<bool> want_primary) {
    if (pid == skip_id) {
      return false;
    }
    Particle *p = cell_structure.get_local_particle(pid);
    if (not p) {
      return false;
    }
    auto &bond_list = p->bonds();
    auto const it =
        find_matching_bond(bond_list, bond_id, sorted_others, want_primary);
    if (it == bond_list.end()) {
      return false;
    }
    bond_list.erase(it);
    return true;
  };

  if (not particle_ids.empty()) {
    auto others0 = other_participant_ids(particle_ids, 0);
    std::ranges::sort(others0);
    // Prefer the unambiguous interpretation: particle_ids[0] holds the
    // primary, every other participant the matching mirror. This
    // disambiguates a bond that exists twice with swapped ownership
    // (added once from each side), where role-independent matching could
    // erase a mismatched pair and leave inconsistent orphans behind.
    if (try_remove(particle_ids[0], others0, true)) {
      for (std::size_t i = 1; i < particle_ids.size(); ++i) {
        auto others = other_participant_ids(particle_ids, i);
        std::ranges::sort(others);
        try_remove(particle_ids[i], others, false);
      }
      return true;
    }

    // Fall back to role-independent matching, e.g. when particle_ids[0]
    // only holds a mirror. Reuse others0 instead of recomputing it.
    bool removed = try_remove(particle_ids[0], others0, std::nullopt);
    for (std::size_t i = 1; i < particle_ids.size(); ++i) {
      auto others = other_participant_ids(particle_ids, i);
      std::ranges::sort(others);
      if (try_remove(particle_ids[i], others, std::nullopt)) {
        removed = true;
      }
    }
    return removed;
  }

  return false;
}

void sync_bond_tuples(std::vector<std::pair<int, std::vector<int>>> &tuples,
                      boost::mpi::communicator const &comm) {
  if (comm.size() > 1) {
    Utils::Mpi::gather_buffer(tuples, comm);
    boost::mpi::broadcast(comm, tuples, 0);
  }
}

void rebuild_bond_mirrors(System::System &system) {
  auto &cell_structure = *system.cell_structure;

  // Gather the primary entries of all ranks, so mirrors can be added for
  // participants living on a different rank than the primary.
  std::vector<std::pair<int, std::vector<int>>> primaries;
  for (auto const &p : cell_structure.local_particles()) {
    for (auto const bond : p.bonds()) {
      if (bond.is_primary()) {
        std::vector<int> ids = {p.id()};
        std::ranges::copy(bond.partner_ids(), std::back_inserter(ids));
        primaries.emplace_back(bond.bond_id(), std::move(ids));
      }
    }
  }
  ::sync_bond_tuples(primaries, ::comm_cart);

  for (auto const &[bond_id, ids] : primaries) {
    // Index 0 is the owner, whose primary entry already exists; only the
    // other participants may need a mirror.
    for (std::size_t i = 1; i < ids.size(); ++i) {
      Particle *p = cell_structure.get_local_particle(ids[i]);
      if (not p) {
        continue;
      }
      auto const others = other_participant_ids(ids, i);
      auto sorted_others = others;
      std::ranges::sort(sorted_others);
      auto const exists = find_matching_bond(p->bonds(), bond_id, sorted_others,
                                             std::nullopt) != p->bonds().end();
      if (not exists) {
        p->bonds().insert(
            BondView(bond_id, {others.data(), others.size()}, false));
      }
    }
  }
}
