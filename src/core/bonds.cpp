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
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/mpi/all_compare.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <utility>
#include <vector>

void add_bond(System::System &system, int bond_id,
              std::vector<int> const &particle_ids) {
  system.bond_storage->add(bond_id, particle_ids);
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
