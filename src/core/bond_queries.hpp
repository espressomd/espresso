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
 * Storage-agnostic, read-only bond queries (see #5046).
 *
 * Code outside of the bond storage implementation (@ref bonds.cpp, particle
 * serialization and ghost communication) should inspect bonds through the
 * helpers taking a @ref Particle, rather than accessing
 * @ref Particle::bonds() and filtering on the storage role of the entries
 * itself, so that the storage can move out of the particle without touching
 * the callers. Mutations go through @ref add_bond() and @ref remove_bond().
 *
 * Everything here is header-only and inlined, since some of it is called
 * from force kernels.
 */

#include "BondList.hpp"
#include "Particle.hpp"

#include <algorithm>
#include <optional>
#include <utility>

/** @brief Which bond entries of a particle to visit. */
enum class BondEntries {
  /**
   * Only the bonds this particle owns. Over all particles, every bond is
   * visited exactly once, with its partner ids in creation order (which
   * angle and dihedral bonds rely on). Use for anything that evaluates,
   * counts or exports bonds.
   */
  owned,
  /**
   * Every bond this particle participates in, whether it owns it or not.
   * A bond is visited once per participant. Use for questions about a
   * particle's connectivity, e.g. which bonds to remove with it.
   */
  all,
};

/**
 * @brief Execute a kernel for the entries of a bond list.
 *
 * Storage-level overload; prefer the overload taking a @ref Particle.
 *
 * @tparam entries Which entries to visit, see @ref BondEntries.
 * @param bonds    Bond list of the particle.
 * @param kernel   Callable invoked as <tt>kernel(BondView const &)</tt>.
 *                 With @ref BondEntries::all, @ref BondView::is_primary()
 *                 tells whether the particle owns the bond.
 */
template <BondEntries entries = BondEntries::owned, typename Kernel>
inline void for_each_bond(BondList const &bonds, Kernel &&kernel) {
  for (auto const bond : bonds) {
    if constexpr (entries == BondEntries::owned) {
      if (not bond.is_primary()) {
        continue;
      }
    }
    kernel(bond);
  }
}

/**
 * @brief Execute a kernel for the bonds of a particle.
 *
 * @tparam entries Which entries to visit, see @ref BondEntries.
 * @param p        Particle whose bonds are visited.
 * @param kernel   Callable invoked as <tt>kernel(BondView const &)</tt>.
 *                 With @ref BondEntries::all, @ref BondView::is_primary()
 *                 tells whether the particle owns the bond.
 */
template <BondEntries entries = BondEntries::owned, typename Kernel>
inline void for_each_bond(Particle const &p, Kernel &&kernel) {
  for_each_bond<entries>(p.bonds(), std::forward<Kernel>(kernel));
}

/**
 * @brief Find the first bond of a particle whose bond id is accepted by
 * @p pred.
 *
 * @tparam entries Which entries to search, see @ref BondEntries.
 * @param p     Particle whose bonds are searched.
 * @param pred  Callable <tt>bool(int bond_id)</tt>.
 * @return The bond id, if any.
 */
template <BondEntries entries = BondEntries::owned, typename Pred>
inline std::optional<int> find_bond_id(Particle const &p, Pred &&pred) {
  for (auto const bond : p.bonds()) {
    if constexpr (entries == BondEntries::owned) {
      if (not bond.is_primary()) {
        continue;
      }
    }
    if (pred(bond.bond_id())) {
      return bond.bond_id();
    }
  }
  return std::nullopt;
}

/**
 * @brief Check whether a pair bond accepted by @p pred connects two
 * particles.
 *
 * Both particles are searched, since one of them may be a ghost whose copy
 * of a recently added bond has not been communicated yet.
 *
 * @param p1    First particle.
 * @param p2    Second particle.
 * @param pred  Callable <tt>bool(int bond_id)</tt> selecting the bonds of
 *              interest.
 */
template <typename Pred>
inline bool pair_bond_exists(Particle const &p1, Particle const &p2,
                             Pred &&pred) {
  if (&p1 == &p2) {
    return false;
  }
  auto const exists_on = [&pred](Particle const &p, int partner_id) {
    auto const &bonds = p.bonds();
    return std::any_of(bonds.begin(), bonds.end(), [&](BondView const &bond) {
      return bond.partner_ids().size() == 1u and
             bond.partner_ids()[0] == partner_id and pred(bond.bond_id());
    });
  };
  return exists_on(p1, p2.id()) or exists_on(p2, p1.id());
}

/**
 * @brief Check whether a pair bond with id @p bond_id connects two
 * particles.
 */
inline bool pair_bond_exists(Particle const &p1, Particle const &p2,
                             int bond_id) {
  return pair_bond_exists(p1, p2, [bond_id](int id) { return id == bond_id; });
}
