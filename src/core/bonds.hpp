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

#pragma once

#include "BondStorage.hpp"
#include "bond_queries.hpp"
#include "system/System.hpp"

#include <boost/mpi/communicator.hpp>

#include <utility>
#include <vector>

/** @brief A bond id and the ids of the bonded particles, owner first. */
using BondTuple = std::pair<int, std::vector<int>>;

/**
 * @brief Add a bond between particles.
 *
 * Must be called on all MPI ranks, with the same arguments and in the same
 * order relative to other bond mutations, since it registers the bond in
 * the global bond topology (@ref System::System::bond_storage).
 *
 * Also writes the primary entry (used for force/energy calculation) on
 * @p particle_ids[0] and a mirror entry on every other participant that is
 * locally known, so the bond can be found and removed from any participant.
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 *
 * @return whether at least one entry was written locally.
 */
bool add_bond(System::System &system, int bond_id,
              std::vector<int> const &particle_ids);

/**
 * @brief Remove a bond between particles.
 *
 * Must be called on all MPI ranks, with the same arguments and in the same
 * order relative to other bond mutations, since it removes the bond from
 * the global bond topology (@ref System::System::bond_storage).
 *
 * Also removes, from every one of @p particle_ids that is locally known, the
 * bond entry (primary or mirror) matching @p bond_id whose partner ids are
 * exactly the remaining ids of @p particle_ids (order does not matter).
 *
 * @param skip_id If not -1, this participant's entry is left untouched, for
 * use when that particle's bond list is about to be discarded anyway (e.g.
 * particle removal).
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 *
 * @return whether at least one entry was removed locally.
 */
bool remove_bond(System::System &system, int bond_id,
                 std::vector<int> const &particle_ids, int skip_id = -1);

/**
 * @brief Recreate the global bond topology and missing mirror entries from
 * the primary entries.
 *
 * Collective over all MPI ranks. The global bond topology is replaced by
 * the primary entries found on all ranks. For every one of them, ensures a
 * matching mirror entry exists on every other participant; existing mirrors
 * are left untouched. To be called after loading bonds that bypassed
 * @ref add_bond(), e.g. from an mpiio checkpoint.
 */
void rebuild_bond_mirrors(System::System &system);

/**
 * @brief Combine bond-removal tuples collected locally on each rank into
 * the same set on every rank.
 *
 * A tuple is @c (bond_id, particle_ids), suitable for a subsequent
 * @ref remove_bond() call. Collecting @p tuples generally only finds the
 * entries local to the calling rank; since @ref remove_bond() must run on
 * every rank to also reach participants living elsewhere, the tuples found
 * on any one rank need to be gathered and broadcast to all ranks first.
 * A no-op when @p comm has a single rank.
 */
void sync_bond_tuples(std::vector<BondTuple> &tuples,
                      boost::mpi::communicator const &comm);

/**
 * @brief Bonds of the global bond topology a particle participates in.
 *
 * Identical on all MPI ranks, unlike the bond entries of the particle,
 * which are only available where the particle is locally known.
 *
 * @param system  The system.
 * @param pid     Particle id.
 * @param pred    Callable <tt>bool(BondStorage::Bond const &)</tt> selecting
 *                the bonds of interest.
 * @return Tuples suitable for @ref remove_bond().
 */
template <typename Pred>
std::vector<BondTuple> get_bonds_of_particle(System::System const &system,
                                             int pid, Pred &&pred) {
  std::vector<BondTuple> bonds;
  system.bond_storage->for_each_bond_of(
      pid, [&](BondStorage::Bond const &bond) {
        if (pred(bond)) {
          bonds.emplace_back(bond.bond_id,
                             std::vector<int>(bond.particle_ids.begin(),
                                              bond.particle_ids.end()));
        }
      });
  return bonds;
}

/**
 * @brief Check that the global bond topology is identical on all MPI ranks
 * and matches the bonds stored on the particles.
 *
 * Collective over all MPI ranks, with a cost proportional to the total
 * number of bonds; meant for tests and additional checks.
 */
bool is_bond_storage_consistent(System::System const &system);
