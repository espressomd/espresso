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

#include "system/System.hpp"

#include <boost/mpi/communicator.hpp>

#include <utility>
#include <vector>

/**
 * @brief Add a bond between particles.
 *
 * Writes the primary entry (used for force/energy calculation) on
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
 * Removes, from every one of @p particle_ids that is locally known, the
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
 * @brief Recreate missing mirror entries from primary entries.
 *
 * Collective over all MPI ranks. For every primary entry found on any rank,
 * ensures a matching mirror entry exists on every other participant.
 * Existing mirrors are left untouched, so this is safe to call
 * unconditionally, e.g. after loading bonds that bypassed @ref add_bond()
 * (an mpiio checkpoint written before bonds were stored on all
 * participants).
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
void sync_bond_tuples(std::vector<std::pair<int, std::vector<int>>> &tuples,
                      boost::mpi::communicator const &comm);
