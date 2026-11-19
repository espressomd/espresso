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

#include "BondStorage.hpp"
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/mpi/all_compare.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/operations.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

void add_bond(System::System &system, int bond_id,
              std::vector<int> const &particle_ids) {
  system.bond_storage->add(bond_id, particle_ids);
}

void check_particles_exist(System::System const &system,
                           std::span<const int> particle_ids) {
  auto const &cell_structure = *system.cell_structure;
  // a particle exists if it resides on one of the ranks
  std::vector<int> found_local;
  found_local.reserve(particle_ids.size());
  for (auto const pid : particle_ids) {
    auto const *p =
        (pid >= 0) ? cell_structure.get_local_particle(pid) : nullptr;
    found_local.push_back((p != nullptr and not p->is_ghost()) ? 1 : 0);
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

static std::string format_ids(std::span<const int> ids) {
  std::string out = "[";
  for (std::size_t i = 0u; i < ids.size(); ++i) {
    out += ((i == 0u) ? "" : ", ") + std::to_string(ids[i]);
  }
  return out + "]";
}

/** @brief Whether @p particle_ids[0] owns exactly this bond. */
static bool owns_bond(System::System const &system, int bond_id,
                      std::span<const int> particle_ids) {
  return system.bond_storage->any_owned_bond(
      particle_ids[0], [&](BondStorage::Bond const &bond) {
        return bond.bond_id == bond_id and
               std::ranges::equal(bond.particle_ids, particle_ids);
      });
}

static void check_not_empty(std::span<const int> particle_ids) {
  if (particle_ids.empty()) {
    throw std::invalid_argument("At least one particle is required");
  }
}

void check_bond_can_be_added(System::System const &system, int bond_id,
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
  check_particles_exist(system, particle_ids);
  if (owns_bond(system, bond_id, particle_ids)) {
    throw std::runtime_error("Bond " + std::to_string(bond_id) +
                             " with partners " + format_ids(partner_ids) +
                             " already exists on particle " +
                             std::to_string(owner_id));
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

bool remove_bond(System::System &system, int bond_id,
                 std::vector<int> const &particle_ids) {
  return system.bond_storage->remove(bond_id, particle_ids);
}

void sync_bond_tuples(std::vector<BondTuple> &tuples,
                      boost::mpi::communicator const &comm) {
  if (comm.size() > 1) {
    Utils::Mpi::gather_buffer(tuples, comm);
    boost::mpi::broadcast(comm, tuples, 0);
  }
}

bool is_bond_storage_consistent(System::System const &system) {
  std::vector<BondTuple> stored;
  system.bond_storage->for_each_bond([&](BondStorage::Bond const &bond) {
    stored.emplace_back(
        bond.bond_id,
        std::vector<int>(bond.particle_ids.begin(), bond.particle_ids.end()));
  });
  std::ranges::sort(stored);
  return Utils::Mpi::all_compare(::comm_cart, stored);
}
