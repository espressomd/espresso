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
#include <iterator>
#include <map>
#include <optional>
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

void remove_bond(System::System &system, int bond_id,
                 std::span<const int> particle_ids) {
  system.cell_structure->bond_table().remove(bond_id, particle_ids);
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

namespace {
/**
 * @brief Checks of a batch of bonds that need communication, with a fixed
 * number of collectives: all particles exist, no bond appears twice in the
 * batch, and no bond exists already.
 */
void check_batch_can_be_added(System::System const &system,
                              std::span<BondRecord const> bonds) {
  // one collective for the existence of all particles
  std::vector<int> unique_ids;
  for (auto const &bond : bonds) {
    std::ranges::copy(bond.particle_ids(), std::back_inserter(unique_ids));
  }
  std::ranges::sort(unique_ids);
  auto const duplicates = std::ranges::unique(unique_ids);
  unique_ids.erase(duplicates.begin(), duplicates.end());
  check_particles_exist(system, unique_ids);
  // no bond twice in the batch
  std::vector<BondRecord> sorted(bonds.begin(), bonds.end());
  std::ranges::sort(sorted, [](BondRecord const &a, BondRecord const &b) {
    return a.bond_id() != b.bond_id() ? a.bond_id() < b.bond_id()
                                      : std::ranges::lexicographical_compare(
                                            a.particle_ids(), b.particle_ids());
  });
  auto const twice = std::ranges::adjacent_find(
      sorted, [](BondRecord const &a, BondRecord const &b) {
        return a.matches(b.bond_id(), b.particle_ids());
      });
  if (twice != sorted.end()) {
    throw std::runtime_error(
        "Bond " + std::to_string(twice->bond_id()) + " with partners " +
        format_ids(twice->partner_ids()) + " appears twice for particle " +
        std::to_string(twice->owner_id()));
  }
  // one collective for the bonds that already exist: their owners know
  auto const exists_here = [&system](BondRecord const &bond) {
    return local_particle(system, bond.owner_id()) != nullptr and
           system.cell_structure->bond_table().contains(bond.bond_id(),
                                                        bond.particle_ids());
  };
  auto first_existing = static_cast<std::size_t>(
      std::ranges::find_if(bonds, exists_here) - bonds.begin());
  first_existing = boost::mpi::all_reduce(::comm_cart, first_existing,
                                          boost::mpi::minimum<std::size_t>());
  if (first_existing != bonds.size()) {
    auto const &bond = bonds[first_existing];
    throw std::runtime_error(
        already_exists_message(bond.bond_id(), bond.particle_ids()));
  }
}
} // namespace

void check_bonds_can_be_added(System::System const &system, int bond_id,
                              std::span<const int> particle_ids,
                              std::size_t n_particles) {
  if (n_particles == 0u or particle_ids.size() % n_particles != 0u) {
    throw std::invalid_argument(
        "Every bond needs the same number of particles");
  }
  std::vector<BondRecord> bonds;
  for (std::size_t i = 0u; i < particle_ids.size(); i += n_particles) {
    bonds.emplace_back(bond_id, particle_ids.subspan(i, n_particles));
    check_bond_shape(system, bond_id, bonds.back().particle_ids());
  }
  check_batch_can_be_added(system, bonds);
}

void add_encoded_bonds(System::System &system, std::span<const int> encoded) {
  // decode and validate all bonds before adding any
  std::vector<BondRecord> bonds;
  while (not encoded.empty()) {
    auto const [bond, rest] = BondEncoding::read(encoded);
    bonds.push_back(bond);
    encoded = rest;
  }
  check_batch_can_be_added(system, bonds);
  for (auto const &bond : bonds) {
    add_bond(system, bond.bond_id(), bond.particle_ids());
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

void check_bonded_interaction_unused(System::System const &system,
                                     std::optional<int> bond_id) {
  auto const rows = system.cell_structure->bond_table().rows();
  auto const used_here = std::ranges::any_of(rows, [&](BondRow const &row) {
    return not row.is_free() and (not bond_id or row.bond_id == *bond_id);
  });
  auto const used =
      boost::mpi::all_reduce(::comm_cart, used_here, std::logical_or<>());
  if (used) {
    throw std::runtime_error(
        bond_id
            ? "Cannot delete bonded interaction " + std::to_string(*bond_id) +
                  ": bonds between particles use it, remove them first"
            : std::string("Cannot delete the bonded interactions: bonds "
                          "between particles use them, remove them "
                          "first"));
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
          bonds.push_back(BondTuple::from(bond));
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
        bonds.push_back(BondTuple::from(bond));
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
      auto const owner_id = bond.particle_ids.front();
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
  // every bond, once per local particle that lists it
  std::vector<std::pair<BondTuple, int>> listings;
  for (auto const &p : system.cell_structure->local_particles()) {
    for (auto const bond : system.cell_structure->bond_table().view(p.id())) {
      listings.emplace_back(BondTuple::from(bond), p.id());
    }
  }
  Utils::Mpi::gather_buffer(listings, ::comm_cart);
  auto consistent = true;
  if (::comm_cart.rank() == 0) {
    // for every bond, the number of times each particle lists it
    std::map<BondTuple, std::map<int, int>> holders;
    for (auto const &[bond, holder] : listings) {
      ++holders[bond][holder];
    }
    for (auto const &[bond, counts] : holders) {
      auto const &ids = bond.particle_ids;
      // each distinct participant lists it the same number of times
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
