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

/**
 * @file
 * Mutations of the bond topology.
 *
 * Every rank stores each bond one of its local particles participates in
 * once, in its @ref BondTable, separate from the particles. A bond is
 * therefore known on every MPI rank where one of its participants resides,
 * and nowhere else. The bonds are shipped along with the particles when
 * they migrate between ranks. Mutations are therefore collective calls in
 * which every rank applies the part that concerns its local particles; no
 * rank holds the complete topology.
 */

#include "BondTable.hpp"
#include "bond_queries.hpp"
#include "system/System.hpp"

#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>

#include <compare>
#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <utility>
#include <vector>

/**
 * @brief A bond that owns its particle ids, e.g. to communicate it, unlike
 * the @ref BondRecord view.
 */
struct BondTuple {
  /** @brief Key of the bonded interaction. */
  int bond_id;
  /** @brief Participants in creation order, owner first. */
  std::vector<int> particle_ids;

  static BondTuple from(BondRecord const &bond) {
    auto const ids = bond.particle_ids();
    return {bond.bond_id(), {ids.begin(), ids.end()}};
  }
  auto operator<=>(BondTuple const &) const = default;

  template <class Archive> void serialize(Archive &ar, long int /*version*/) {
    ar & bond_id;
    ar & particle_ids;
  }
};

/**
 * @brief Add a bond between particles.
 *
 * Must be called on all MPI ranks with the same arguments. Every rank
 * where a participant resides stores the bond. No communication takes
 * place.
 *
 * The first particle owns the bond: the bond is evaluated on the MPI rank
 * where that particle resides. Bonds given by the user should be validated
 * with @ref check_bond_can_be_added() first.
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 */
void add_bond(System::System &system, int bond_id,
              std::span<const int> particle_ids);

/** @overload */
inline void add_bond(System::System &system, int bond_id,
                     std::vector<int> const &particle_ids) {
  add_bond(system, bond_id, std::span<const int>(particle_ids));
}

/**
 * @brief Check that particles exist, e.g. before bonding them.
 *
 * Collective over all MPI ranks, since a rank only knows the particles in
 * its own domain.
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
 * Collective over all MPI ranks.
 *
 * @throw std::invalid_argument on all ranks, if the bond is malformed.
 * @throw std::runtime_error on all ranks, if a particle does not exist or
 *        the bond already exists.
 */
void check_bond_can_be_added(System::System const &system, int bond_id,
                             std::span<const int> particle_ids);

/**
 * @brief Check that bonds can be added, e.g. before a user adds them.
 *
 * Same checks as @ref check_bond_can_be_added() for every bond, plus no
 * bond may appear twice in the batch. Uses a fixed number of collectives,
 * independent of the number of bonds.
 *
 * @param system        The system.
 * @param bond_id       Key of the bonded interaction.
 * @param particle_ids  Participants of all bonds, one bond after the other,
 *                      each with the owner first.
 * @param n_particles   Number of participants per bond.
 * @throw std::invalid_argument on all ranks, if a bond is malformed.
 * @throw std::runtime_error on all ranks, if a particle does not exist or
 *        a bond already exists.
 */
void check_bonds_can_be_added(System::System const &system, int bond_id,
                              std::span<const int> particle_ids,
                              std::size_t n_particles);

/**
 * @brief Add bonds of the same bonded interaction between particles.
 *
 * Equivalent to calling @ref add_bond() for every bond. Must be called on
 * all MPI ranks with the same arguments; no communication takes place.
 *
 * @param system        The system.
 * @param bond_id       Key of the bonded interaction.
 * @param particle_ids  Participants of all bonds, see
 *                      @ref check_bonds_can_be_added().
 * @param n_particles   Number of participants per bond.
 */
void add_bonds(System::System &system, int bond_id,
               std::span<const int> particle_ids, std::size_t n_particles);

/**
 * @brief Add bonds of any bonded interactions between existing particles.
 *
 * The bonds are encoded one after the other as
 * <tt>[bond_id, n, id_0, ..., id_{n-1}]</tt> and added in this order, as
 * with @ref add_bond(). Only the existence of the particles is checked,
 * with a single collective. Must be called on all MPI ranks with the same
 * arguments.
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 *
 * @throw std::invalid_argument on all ranks, if the encoding is malformed.
 * @throw std::runtime_error on all ranks, if a particle does not exist;
 *        no bond is added then.
 */
void add_encoded_bonds(System::System &system, std::span<const int> encoded);

/**
 * @brief Check that a bond can be removed, e.g. before a user removes it.
 *
 * The bond must be owned by @p particle_ids[0], with the remaining
 * particles as its partners, in this order.
 *
 * Collective over all MPI ranks.
 *
 * @throw std::invalid_argument on all ranks, if @p particle_ids is empty.
 * @throw std::runtime_error on all ranks, if the bond does not exist.
 */
void check_bond_can_be_removed(System::System const &system, int bond_id,
                               std::span<const int> particle_ids);

/**
 * @brief Remove a bond between particles.
 *
 * Must be called on all MPI ranks with the same arguments. Every rank
 * removes the bond matching @p bond_id and exactly @p particle_ids, in the
 * given order, if it stores it. Nothing happens if the bond does not
 * exist; bonds given by the user should be validated with
 * @ref check_bond_can_be_removed() first. No communication takes place.
 *
 * The caller is responsible for calling
 * @ref System::System::on_particle_change().
 */
void remove_bond(System::System &system, int bond_id,
                 std::span<const int> particle_ids);

/** @overload */
inline void remove_bond(System::System &system, int bond_id,
                        std::vector<int> const &particle_ids) {
  remove_bond(system, bond_id, std::span<const int>(particle_ids));
}

/**
 * @brief Check that no bond uses a bonded interaction, e.g. before the
 * bonded interaction is deleted.
 *
 * Collective over all MPI ranks.
 *
 * @param system   The system.
 * @param bond_id  Key of the bonded interaction, or nothing to check that
 *                 there are no bonds at all.
 * @throw std::runtime_error on all ranks, if a bond uses the bonded
 *        interaction.
 */
void check_bonded_interaction_unused(System::System const &system,
                                     std::optional<int> bond_id);

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
 * @brief Collect bonds of particles, identically on all MPI ranks.
 *
 * Each rank collects the bonds of those particles in @p pids that reside
 * on it, then the lists are combined (see @ref sync_bond_tuples()).
 * Collective over all MPI ranks.
 *
 * @param system  The system.
 * @param pids    Particle ids.
 * @param pred    Callable <tt>bool(BondRecord const &, int pid)</tt>
 *                selecting the bonds of interest of particle @c pid.
 *                A bond is collected once per particle that selects it.
 * @return Tuples suitable for @ref remove_bond().
 */
std::vector<BondTuple> gather_bonds_of_particles(
    System::System const &system, std::span<const int> pids,
    std::function<bool(BondRecord const &, int)> const &pred);

/**
 * @brief Fetch the bonds owned by particles from the ranks where they
 * reside, with a single collective, into a cache on the head node.
 *
 * Head node only, for code that runs on the head node alone, e.g. the
 * analysis of particle copies fetched with @ref get_particle_data(). The
 * cache is cleared by @ref invalidate_bond_fetch_cache(), which is called
 * whenever particles change.
 */
void prefetch_owned_bonds(std::vector<int> const &pids);

/**
 * @brief Bonds owned by a particle, fetched from the rank where it resides.
 *
 * Head node only, see @ref prefetch_owned_bonds(). The reference is valid
 * until the next fetch or invalidation.
 */
std::vector<BondTuple> const &get_owned_bonds(int pid);

/** @brief Clear the cache of @ref get_owned_bonds(). */
void invalidate_bond_fetch_cache();

/**
 * @brief Check that the bonds are listed consistently.
 *
 * Every bond must be listed by each of its distinct participants, as seen
 * from the rank where the participant resides, the same number of times,
 * and by no other particle. Collective over all MPI
 * ranks, with a cost proportional to the total number of bonds; meant for
 * tests and additional checks.
 *
 * @return the result, on all ranks.
 */
bool is_bond_storage_consistent(System::System const &system);
