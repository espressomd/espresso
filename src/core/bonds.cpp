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

#include "bonds.hpp"

#include "BondTable.hpp"
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/operations.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
/** @brief The particle with id @p pid, if it resides on this rank. */
Particle *local_particle(System::System const &system, int pid) {
  if (pid < 0) {
    return nullptr;
  }
  auto *p = system.cell_structure->get_local_particle(pid);
  return (p != nullptr and not p->is_ghost()) ? p : nullptr;
}

/** @brief The distinct participants of a bond residing on this rank. */
std::vector<int> local_participants(System::System const &system,
                                    std::span<const int> particle_ids) {
  std::vector<int> pids;
  for (auto it = particle_ids.begin(); it != particle_ids.end(); ++it) {
    if (std::find(particle_ids.begin(), it, *it) == it and
        local_particle(system, *it) != nullptr) {
      pids.push_back(*it);
    }
  }
  return pids;
}

std::string format_ids(std::span<const int> ids) {
  std::string out = "[";
  for (std::size_t i = 0u; i < ids.size(); ++i) {
    out += ((i == 0u) ? "" : ", ") + std::to_string(ids[i]);
  }
  return out + "]";
}

/**
 * @brief Whether @p particle_ids[0] owns exactly this bond.
 * Collective: only the rank where the owner resides knows the answer.
 */
bool owns_bond(System::System const &system, int bond_id,
               std::span<const int> particle_ids) {
  auto const local =
      local_particle(system, particle_ids[0]) != nullptr and
      system.cell_structure->bond_table().contains(bond_id, particle_ids);
  return boost::mpi::all_reduce(::comm_cart, local, std::logical_or<>());
}

void check_not_empty(std::span<const int> particle_ids) {
  if (particle_ids.empty()) {
    throw std::invalid_argument("At least one particle is required");
  }
}
} // namespace

void add_bond(System::System &system, int bond_id,
              std::span<const int> particle_ids) {
  system.cell_structure->bond_table().add(
      bond_id, particle_ids, local_participants(system, particle_ids));
}

bool remove_bond(System::System &system, int bond_id,
                 std::span<const int> particle_ids) {
  if (particle_ids.empty()) {
    return false;
  }
  auto const removed =
      system.cell_structure->bond_table().remove(bond_id, particle_ids);
  return removed and local_particle(system, particle_ids[0]) != nullptr;
}

void check_particles_exist(System::System const &system,
                           std::span<const int> particle_ids) {
  // a particle exists if it resides on one of the ranks
  std::vector<int> found_local;
  found_local.reserve(particle_ids.size());
  for (auto const pid : particle_ids) {
    found_local.push_back((local_particle(system, pid) != nullptr) ? 1 : 0);
  }
  auto found = found_local;
  if (not particle_ids.empty()) {
    boost::mpi::all_reduce(::comm_cart, found_local.data(),
                           static_cast<int>(found_local.size()), found.data(),
                           boost::mpi::maximum<int>());
  }
  for (std::size_t i = 0u; i < particle_ids.size(); ++i) {
    if (found[i] == 0) {
      throw std::runtime_error("Particle with id " +
                               std::to_string(particle_ids[i]) +
                               " does not exist");
    }
  }
}

namespace {
/** @brief Checks of a bond that need no communication. */
void check_bond_shape(System::System const &system, int bond_id,
                      std::span<const int> particle_ids) {
  check_not_empty(particle_ids);
  auto const owner_id = particle_ids[0];
  auto const partner_ids = particle_ids.subspan(1u);
  if (std::ranges::find(partner_ids, owner_id) != partner_ids.end()) {
    throw std::invalid_argument("Bond partners " + format_ids(partner_ids) +
                                " include the particle " +
                                std::to_string(owner_id) + " itself");
  }
  std::vector<int> sorted(partner_ids.begin(), partner_ids.end());
  std::ranges::sort(sorted);
  if (std::ranges::adjacent_find(sorted) != sorted.end()) {
    throw std::invalid_argument("Cannot add duplicate bond partners " +
                                format_ids(partner_ids) + " to particle " +
                                std::to_string(owner_id));
  }
  auto const n_partners = number_of_partners(*system.bonded_ias->at(bond_id));
  if (std::cmp_not_equal(partner_ids.size(), n_partners)) {
    throw std::invalid_argument("Bond " + std::to_string(bond_id) + " needs " +
                                std::to_string(n_partners) + " partners, got " +
                                std::to_string(partner_ids.size()));
  }
}

std::string already_exists_message(int bond_id,
                                   std::span<const int> particle_ids) {
  return "Bond " + std::to_string(bond_id) + " with partners " +
         format_ids(particle_ids.subspan(1u)) + " already exists on particle " +
         std::to_string(particle_ids[0]);
}
} // namespace

void check_bond_can_be_added(System::System const &system, int bond_id,
                             std::span<const int> particle_ids) {
  check_bond_shape(system, bond_id, particle_ids);
  check_particles_exist(system, particle_ids);
  if (owns_bond(system, bond_id, particle_ids)) {
    throw std::runtime_error(already_exists_message(bond_id, particle_ids));
  }
}

void check_bonds_can_be_added(System::System const &system, int bond_id,
                              std::span<const int> particle_ids,
                              std::size_t n_particles) {
  if (n_particles == 0u or particle_ids.size() % n_particles != 0u) {
    throw std::invalid_argument(
        "Every bond needs the same number of particles");
  }
  auto const n_bonds = particle_ids.size() / n_particles;
  auto const bond = [&](std::size_t i) {
    return particle_ids.subspan(i * n_particles, n_particles);
  };
  for (std::size_t i = 0u; i < n_bonds; ++i) {
    check_bond_shape(system, bond_id, bond(i));
  }
  // one collective for the existence of all particles
  std::vector<int> unique_ids(particle_ids.begin(), particle_ids.end());
  std::ranges::sort(unique_ids);
  auto const duplicates = std::ranges::unique(unique_ids);
  unique_ids.erase(duplicates.begin(), duplicates.end());
  check_particles_exist(system, unique_ids);
  // no bond twice in the batch
  std::vector<std::size_t> order(n_bonds);
  for (std::size_t i = 0u; i < n_bonds; ++i) {
    order[i] = i;
  }
  auto const less = [&](std::size_t a, std::size_t b) {
    return std::ranges::lexicographical_compare(bond(a), bond(b));
  };
  std::ranges::sort(order, less);
  for (std::size_t i = 1u; i < n_bonds; ++i) {
    if (std::ranges::equal(bond(order[i - 1u]), bond(order[i]))) {
      throw std::runtime_error(
          "Bond " + std::to_string(bond_id) + " with partners " +
          format_ids(bond(order[i]).subspan(1u)) +
          " appears twice for particle " + std::to_string(bond(order[i])[0]));
    }
  }
  // one collective for the bonds that already exist: their owners know
  auto first_existing = n_bonds;
  for (std::size_t i = 0u; i < n_bonds; ++i) {
    auto const ids = bond(i);
    if (local_particle(system, ids[0]) != nullptr and
        system.cell_structure->bond_table().contains(bond_id, ids)) {
      first_existing = i;
      break;
    }
  }
  first_existing = boost::mpi::all_reduce(::comm_cart, first_existing,
                                          boost::mpi::minimum<std::size_t>());
  if (first_existing != n_bonds) {
    throw std::runtime_error(
        already_exists_message(bond_id, bond(first_existing)));
  }
}

void add_encoded_bonds(System::System &system, std::span<const int> encoded) {
  // decode and validate the layout before touching any particle
  std::vector<std::pair<int, std::span<const int>>> bonds;
  std::vector<int> unique_ids;
  for (std::size_t i = 0u; i < encoded.size();) {
    if (i + 2u > encoded.size() or encoded[i + 1u] <= 0 or
        static_cast<std::size_t>(encoded[i + 1u]) > BondTable::max_particles or
        i + 2u + static_cast<std::size_t>(encoded[i + 1u]) > encoded.size()) {
      throw std::invalid_argument("Malformed bond encoding");
    }
    auto const n = static_cast<std::size_t>(encoded[i + 1u]);
    auto const ids = encoded.subspan(i + 2u, n);
    bonds.emplace_back(encoded[i], ids);
    unique_ids.insert(unique_ids.end(), ids.begin(), ids.end());
    i += 2u + n;
  }
  std::ranges::sort(unique_ids);
  auto const duplicates = std::ranges::unique(unique_ids);
  unique_ids.erase(duplicates.begin(), duplicates.end());
  check_particles_exist(system, unique_ids);
  // a bond can be stored only once
  auto sorted = bonds;
  std::ranges::sort(sorted, [](auto const &a, auto const &b) {
    return a.first != b.first
               ? a.first < b.first
               : std::ranges::lexicographical_compare(a.second, b.second);
  });
  for (std::size_t i = 1u; i < sorted.size(); ++i) {
    if (sorted[i - 1u].first == sorted[i].first and
        std::ranges::equal(sorted[i - 1u].second, sorted[i].second)) {
      throw std::runtime_error(
          "Bond " + std::to_string(sorted[i].first) + " with partners " +
          format_ids(sorted[i].second.subspan(1u)) +
          " appears twice for particle " + std::to_string(sorted[i].second[0]));
    }
  }
  auto first_existing = bonds.size();
  for (std::size_t i = 0u; i < bonds.size(); ++i) {
    auto const &[bond_id, ids] = bonds[i];
    if (local_particle(system, ids[0]) != nullptr and
        system.cell_structure->bond_table().contains(bond_id, ids)) {
      first_existing = i;
      break;
    }
  }
  first_existing = boost::mpi::all_reduce(::comm_cart, first_existing,
                                          boost::mpi::minimum<std::size_t>());
  if (first_existing != bonds.size()) {
    throw std::runtime_error(already_exists_message(
        bonds[first_existing].first, bonds[first_existing].second));
  }
  for (auto const &[bond_id, ids] : bonds) {
    add_bond(system, bond_id, ids);
  }
}

void add_bonds(System::System &system, int bond_id,
               std::span<const int> particle_ids, std::size_t n_particles) {
  assert(n_particles != 0u and particle_ids.size() % n_particles == 0u);
  for (std::size_t i = 0u; i < particle_ids.size(); i += n_particles) {
    add_bond(system, bond_id, particle_ids.subspan(i, n_particles));
  }
}

void check_bond_can_be_removed(System::System const &system, int bond_id,
                               std::span<const int> particle_ids) {
  check_not_empty(particle_ids);
  if (not owns_bond(system, bond_id, particle_ids)) {
    throw std::runtime_error(
        "Bond " + std::to_string(bond_id) + " with partners " +
        format_ids(particle_ids.subspan(1u)) + " doesn't exist on particle " +
        std::to_string(particle_ids[0]));
  }
}

void sync_bond_tuples(std::vector<BondTuple> &tuples,
                      boost::mpi::communicator const &comm) {
  if (comm.size() > 1) {
    Utils::Mpi::gather_buffer(tuples, comm);
    boost::mpi::broadcast(comm, tuples, 0);
  }
}

std::vector<BondTuple> gather_bonds_of_particles(
    System::System const &system, std::span<const int> pids,
    std::function<bool(BondRecord const &, int)> const &pred) {
  std::vector<BondTuple> bonds;
  for (auto const pid : pids) {
    if (local_particle(system, pid) != nullptr) {
      for (auto const bond : system.cell_structure->bond_table().view(pid)) {
        if (pred(bond, pid)) {
          auto const ids = bond.particle_ids();
          bonds.emplace_back(bond.bond_id(),
                             std::vector<int>(ids.begin(), ids.end()));
        }
      }
    }
  }
  sync_bond_tuples(bonds, ::comm_cart);
  return bonds;
}

namespace {
/** @brief Owned bonds of particles, fetched to the head node. */
std::unordered_map<int, std::vector<BondTuple>> bond_fetch_cache;
/** @brief Number of particles above which the cache is cleared. */
constexpr std::size_t bond_fetch_cache_max_size = std::size_t{1} << 20;
} // namespace

static void mpi_gather_owned_bonds_local(std::vector<int> const &pids) {
  auto const &system = System::get_system();
  auto const &bond_table = system.cell_structure->bond_table();
  std::vector<BondTuple> bonds;
  for (auto const pid : pids) {
    if (local_particle(system, pid) != nullptr) {
      for (auto const bond : bond_table.owned(pid)) {
        auto const ids = bond.particle_ids();
        bonds.emplace_back(bond.bond_id(),
                           std::vector<int>(ids.begin(), ids.end()));
      }
    }
  }
  Utils::Mpi::gather_buffer(bonds, ::comm_cart);
  if (::comm_cart.rank() == 0) {
    if (bond_fetch_cache.size() + pids.size() > bond_fetch_cache_max_size) {
      bond_fetch_cache.clear();
    }
    for (auto const pid : pids) {
      bond_fetch_cache[pid].clear();
    }
    for (auto &bond : bonds) {
      auto const owner_id = bond.second.front();
      bond_fetch_cache[owner_id].emplace_back(std::move(bond));
    }
  }
}

REGISTER_CALLBACK(mpi_gather_owned_bonds_local)

void prefetch_owned_bonds(std::vector<int> const &pids) {
  std::vector<int> missing;
  for (auto const pid : pids) {
    if (not bond_fetch_cache.contains(pid)) {
      missing.push_back(pid);
    }
  }
  if (not missing.empty()) {
    Communication::mpiCallbacks().call_all(mpi_gather_owned_bonds_local,
                                           missing);
  }
}

std::vector<BondTuple> const &get_owned_bonds(int pid) {
  prefetch_owned_bonds({pid});
  return bond_fetch_cache.at(pid);
}

void invalidate_bond_fetch_cache() { bond_fetch_cache.clear(); }

bool is_bond_storage_consistent(System::System const &system) {
  // every stored copy, with the id of the particle that holds it
  std::vector<std::pair<BondTuple, int>> copies;
  for (auto const &p : system.cell_structure->local_particles()) {
    for (auto const bond : system.cell_structure->bond_table().view(p.id())) {
      auto const ids = bond.particle_ids();
      copies.emplace_back(
          BondTuple{bond.bond_id(), std::vector<int>(ids.begin(), ids.end())},
          p.id());
    }
  }
  Utils::Mpi::gather_buffer(copies, ::comm_cart);
  auto consistent = true;
  if (::comm_cart.rank() == 0) {
    // for every bond, the number of copies held by each particle
    std::map<BondTuple, std::map<int, int>> holders;
    for (auto const &[bond, holder] : copies) {
      ++holders[bond][holder];
    }
    for (auto const &[bond, counts] : holders) {
      auto const &ids = bond.second;
      // each distinct participant holds the same number of copies
      std::map<int, int> expected;
      for (auto const pid : ids) {
        expected[pid] = counts.begin()->second;
      }
      if (counts != expected) {
        consistent = false;
        break;
      }
    }
  }
  boost::mpi::broadcast(::comm_cart, consistent, 0);
  return consistent;
}
