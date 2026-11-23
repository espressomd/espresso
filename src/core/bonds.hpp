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

#include <span>
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
 * The first particle owns the bond: the bond is evaluated on the MPI rank
 * where that particle resides. The particles need not be known on the
 * calling rank, but bonds given by the user should be validated with
 * @ref check_particles_exist() first.
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 */
void add_bond(System::System &system, int bond_id,
              std::vector<int> const &particle_ids);

/**
 * @brief Check that particles exist, e.g. before bonding them.
 *
 * Collective over all MPI ranks, since a rank only knows the particles in
 * its own domain. @ref add_bond() does not carry out this check itself: it
 * is also called from the integration loop, with particles that are known
 * to exist.
 *
 * @param system        The system.
 * @param particle_ids  Ids of the particles to look for.
 * @throw std::runtime_error on all ranks, if one of the particles does not
 *        exist on any rank.
 */
void check_particles_exist(System::System const &system,
                           std::span<const int> particle_ids);

/**
 * @brief Check that a bond can be added, e.g. before a user adds it.
 *
 * The bond must have the number of partners its bonded interaction
 * expects, its partners must be distinct and differ from the owner
 * (@p particle_ids[0]), all particles must exist, and the owner must not
 * own the same bond yet, with the same partners in the same order.
 *
 * Collective over all MPI ranks, see @ref check_particles_exist().
 *
 * @throw std::invalid_argument on all ranks, if the bond is malformed.
 * @throw std::runtime_error on all ranks, if a particle does not exist or
 *        the bond already exists.
 */
void check_bond_can_be_added(System::System const &system, int bond_id,
                             std::span<const int> particle_ids);

/**
 * @brief Check that a bond can be removed, e.g. before a user removes it.
 *
 * The bond must be owned by @p particle_ids[0], with the remaining
 * particles as its partners, in this order.
 *
 * @throw std::invalid_argument on all ranks, if @p particle_ids is empty.
 * @throw std::runtime_error on all ranks, if the bond does not exist.
 */
void check_bond_can_be_removed(System::System const &system, int bond_id,
                               std::span<const int> particle_ids);

/**
 * @brief Remove a bond between particles.
 *
 * Must be called on all MPI ranks, with the same arguments and in the same
 * order relative to other bond mutations, since it removes the bond from
 * the global bond topology (@ref System::System::bond_storage).
 *
 * Removes one bond matching @p bond_id that is owned by @p particle_ids[0]
 * and whose partners are exactly the remaining ids, in the given order.
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 *
 * @return whether a bond was removed.
 */
bool remove_bond(System::System &system, int bond_id,
                 std::vector<int> const &particle_ids);

/**
 * @brief Combine bond tuples collected locally on each rank into the same
 * list on every rank.
 *
 * A tuple is @c (bond_id, particle_ids), suitable for a subsequent
 * @ref add_bond() or @ref remove_bond() call. Since these must be called
 * identically on every rank, tuples that are only known to some ranks
 * need to be gathered and broadcast to all ranks first.
 * A no-op when @p comm has a single rank.
 */
void sync_bond_tuples(std::vector<BondTuple> &tuples,
                      boost::mpi::communicator const &comm);

/**
 * @brief Bonds of the global bond topology a particle participates in.
 *
 * Identical on all MPI ranks.
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
  system.bond_storage->for_each_bond_with(
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
 * @brief Check that the global bond topology is identical on all MPI ranks.
 *
 * Collective over all MPI ranks, with a cost proportional to the total
 * number of bonds; meant for tests and additional checks.
 */
bool is_bond_storage_consistent(System::System const &system);
