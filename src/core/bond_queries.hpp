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

#pragma once

/**
 * @file
 * Read-only queries on the bonds of a particle.
 *
 * The bonds are read from the global bond topology (@ref BondStorage),
 * which is identical on all MPI ranks: the queries give the same answer
 * everywhere, whether or not the particles are known on the calling rank.
 * Mutations go through @ref add_bond() and @ref remove_bond().
 *
 * Everything here is header-only and inlined, since some of it is called
 * from force kernels. The queries are safe to call concurrently.
 */

#include "BondStorage.hpp"

#include <optional>
#include <utility>

/**
 * @brief Find the first bond owned by a particle whose bond id is accepted
 * by @p pred.
 *
 * @param bond_storage  Global bond topology.
 * @param pid           Particle id.
 * @param pred          Callable <tt>bool(int bond_id)</tt>.
 * @return The bond id, if any.
 */
template <typename Pred>
inline std::optional<int> find_owned_bond_id(BondStorage const &bond_storage,
                                             int pid, Pred &&pred) {
  std::optional<int> result;
  bond_storage.any_owned_bond(pid, [&](BondStorage::Bond const &bond) {
    if (pred(bond.bond_id)) {
      result = bond.bond_id;
    }
    return result.has_value();
  });
  return result;
}

/**
 * @brief Find the first bond a particle participates in, whether it owns
 * it or not, whose bond id is accepted by @p pred.
 *
 * @param bond_storage  Global bond topology.
 * @param pid           Particle id.
 * @param pred          Callable <tt>bool(int bond_id)</tt>.
 * @return The bond id, if any.
 */
template <typename Pred>
inline std::optional<int> find_bond_id_with(BondStorage const &bond_storage,
                                            int pid, Pred &&pred) {
  std::optional<int> result;
  bond_storage.any_bond_with(pid, [&](BondStorage::Bond const &bond) {
    if (pred(bond.bond_id)) {
      result = bond.bond_id;
    }
    return result.has_value();
  });
  return result;
}

/**
 * @brief Check whether a pair bond accepted by @p pred connects two
 * particles, whichever of them owns it.
 *
 * @param bond_storage  Global bond topology.
 * @param pid1          Id of the first particle.
 * @param pid2          Id of the second particle.
 * @param pred          Callable <tt>bool(int bond_id)</tt> selecting the
 *                      bonds of interest.
 */
template <typename Pred>
inline bool pair_bond_exists(BondStorage const &bond_storage, int pid1,
                             int pid2, Pred &&pred) {
  if (pid1 == pid2) {
    return false;
  }
  // A pair bond is owned by one of its two particles: look at the bonds
  // each of them owns, not at all bonds they participate in.
  auto const owns_pair_with = [&](int owner_id, int partner_id) {
    return bond_storage.any_owned_bond(
        owner_id, [&](BondStorage::Bond const &bond) {
          return bond.particle_ids.size() == 2u and
                 bond.particle_ids[1] == partner_id and pred(bond.bond_id);
        });
  };
  return owns_pair_with(pid1, pid2) or owns_pair_with(pid2, pid1);
}

/**
 * @brief Check whether a pair bond with id @p bond_id connects two
 * particles.
 */
inline bool pair_bond_exists(BondStorage const &bond_storage, int pid1,
                             int pid2, int bond_id) {
  return pair_bond_exists(bond_storage, pid1, pid2,
                          [bond_id](int id) { return id == bond_id; });
}
