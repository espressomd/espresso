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

#include <vector>

/**
 * @brief Add a bond between particles.
 *
 * Writes a primary bond entry on @p particle_ids[0] (the entry used for
 * force/energy calculation, holding the remaining ids in the given order)
 * and a mirror entry on every other participant that is locally known
 * (holding the remaining ids in the same relative order, minus itself), so
 * the bond can be found, queried and removed starting from any participant.
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
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 *
 * @return whether at least one entry was removed locally.
 */
bool remove_bond(System::System &system, int bond_id,
                 std::vector<int> const &particle_ids);

/**
 * @brief Recreate missing mirror entries from primary entries.
 *
 * Collective over all MPI ranks. For every primary bond entry found on any
 * rank, ensures a matching mirror entry exists on every other participant
 * that is locally known anywhere in the simulation. Existing mirror entries
 * are left untouched, so this is safe to call unconditionally (e.g. after
 * loading particle data whose bonds were populated directly rather than
 * through @ref add_bond(), such as an mpiio checkpoint written before
 * bonds were stored on all participants).
 */
void rebuild_bond_mirrors(System::System &system);
