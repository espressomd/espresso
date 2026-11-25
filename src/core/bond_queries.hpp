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
 * The bonds of the particles residing on an MPI rank are kept in the
 * @ref BondTable of that rank, separate from the particles: it holds every
 * bond with a participant on the rank, whether that particle owns it or
 * not. Ghost particles have no entry. Mutations go through @ref add_bond() and
 * @ref remove_bond().
 *
 * Everything here is header-only and inlined, since some of it is called
 * from force kernels. The queries are safe to call concurrently.
 */

#include "BondTable.hpp"
#include "Particle.hpp"

#include <optional>
#include <utility>

/**
 * @brief Find the first bond owned by a particle whose bond id is accepted
 * by @p pred.
 *
 * @param bond_table  Bonds of the local particles.
 * @param p           Local particle.
 * @param pred        Callable <tt>bool(int bond_id)</tt>.
 * @return The bond id, if any.
 */
template <typename Pred>
inline std::optional<int> find_owned_bond_id(BondTable const &bond_table,
                                             Particle const &p, Pred &&pred) {
  for (auto const bond : bond_table.owned(p.id())) {
    if (pred(bond.bond_id())) {
      return bond.bond_id();
    }
  }
  return std::nullopt;
}

/**
 * @brief Find the first bond a particle participates in, whether it owns
 * it or not, whose bond id is accepted by @p pred.
 *
 * @param bond_table  Bonds of the local particles.
 * @param p           Local particle.
 * @param pred        Callable <tt>bool(int bond_id)</tt>.
 * @return The bond id, if any.
 */
template <typename Pred>
inline std::optional<int> find_bond_id_with(BondTable const &bond_table,
                                            Particle const &p, Pred &&pred) {
  for (auto const bond : bond_table.view(p.id())) {
    if (pred(bond.bond_id())) {
      return bond.bond_id();
    }
  }
  return std::nullopt;
}

/**
 * @brief Check whether a pair bond accepted by @p pred connects two
 * particles, given the bonds of one of them.
 *
 * @param bonds     All bonds of the particle with id @p pid1, as a range
 *                  of @ref BondRecord.
 * @param pid1      Id of the particle the bonds belong to.
 * @param pid2      Id of the other particle.
 * @param pred      Callable <tt>bool(int bond_id)</tt>.
 */
template <typename Bonds, typename Pred>
inline bool pair_bond_exists(Bonds const &bonds, int pid1, int pid2,
                             Pred &&pred) {
  if (pid1 == pid2) {
    return false;
  }
  for (auto const bond : bonds) {
    auto const ids = bond.particle_ids();
    if (ids.size() == 2u and (ids[0] == pid2 or ids[1] == pid2) and
        pred(bond.bond_id())) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Check whether a pair bond accepted by @p pred connects two
 * particles, whichever of them owns it.
 *
 * The bonds are read from the first particle that is not a ghost, since
 * ghosts have no entry in the bond table. In pair loops, at least one of
 * the two particles is local.
 *
 * @param bond_table  Bonds of the local particles.
 * @param p1          First particle.
 * @param p2          Second particle.
 * @param pred        Callable <tt>bool(int bond_id)</tt> selecting the bonds
 *                    of interest.
 */
template <typename Pred>
inline bool pair_bond_exists(BondTable const &bond_table, Particle const &p1,
                             Particle const &p2, Pred &&pred) {
  auto const &holder = (p1.is_ghost()) ? p2 : p1;
  auto const &other = (p1.is_ghost()) ? p1 : p2;
  return pair_bond_exists(bond_table.view(holder.id()), holder.id(), other.id(),
                          std::forward<Pred>(pred));
}

/**
 * @brief Check whether a pair bond with id @p bond_id connects two
 * particles.
 */
inline bool pair_bond_exists(BondTable const &bond_table, Particle const &p1,
                             Particle const &p2, int bond_id) {
  return pair_bond_exists(bond_table, p1, p2,
                          [bond_id](int id) { return id == bond_id; });
}
