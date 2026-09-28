/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "cell_system/AtomDecomposition.hpp"
#include "cell_system/ParticleDecomposition.hpp"
#include "cell_system/RegularDecomposition.hpp"

#include "cell_system/Cell.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <cstddef>
#include <functional>
#include <optional>
#include <set>
#include <span>
#include <utility>
#include <vector>

/**
 * @brief Hybrid decomposition cell system.
 *
 * Store particles with short-range interactions
 * in a @ref RegularDecomposition cell system and
 * particles with long-range interactions
 * in a @ref AtomDecomposition (N-square) cell system.
 * All regular cells are coupled to the N-square cells.
 */
class HybridDecomposition : public ParticleDecomposition {
  boost::mpi::communicator m_comm;
  BoxGeometry const &m_box;
  double m_cutoff_regular;
  std::vector<Cell *> m_local_cells;
  std::vector<Cell *> m_ghost_cells;

  /**
   * Topology-agnostic direct-neighbor halo plan (see @c make_halo_plan).
   * Holds ParticleList pointers into this decomposition's cells.
   * Value-copying this object leaves these pointers dangling.
   * @todo make non-copyable or rebuild-on-copy.
   */
  GhostComm::HaloPlan m_halo_plan;

  /** RegularDecomposition to hold the small particles */
  RegularDecomposition m_regular_decomposition;
  /** N-Square Decomposition to hold large particles */
  AtomDecomposition m_n_square;
  /** Set containing the types that should be handled using n_square */
  std::set<int> const m_n_square_types;

  std::function<bool()> m_get_global_ghost_flags;

  /** Commit the ParticleStore: child resorts stage migrated particles into
   *  cells without committing them to store rows; the internal ghost
   *  communications in @ref resort need the committed rows/columns. This
   *  callback (set by @ref CellStructure via @ref set_commit_store) runs
   *  @ref CellStructure::ensure_particle_store_synchronized. */
  std::function<void()> m_commit_store;

  /** Shared migration staging handle: the hybrid's own type-based moves between
   *  the two child decompositions copy a live row into this staging store and
   *  stage a reference to it, exactly like a child's wrong-cell-local move; the
   *  children hold the SAME handle (same address-stable staging-store member on
   *  @ref CellStructure). */
  MigrationStaging m_migration_staging;

  bool is_n_square_type(int type_id) const {
    return m_n_square_types.contains(type_id);
  }

public:
  HybridDecomposition(boost::mpi::communicator comm, double cutoff_regular,
                      double skin, std::function<bool()> get_ghost_flags,
                      BoxGeometry const &box_geo, LocalBox const &local_box,
                      std::set<int> n_square_types);

  /** @brief Set the store-commit callback; see @c m_commit_store. */
  void set_commit_store(std::function<void()> commit) {
    m_commit_store = std::move(commit);
  }

  /** @brief Install the migration staging handle and propagate it to the two
   *  child decompositions, which run the actual wire exchange in their
   *  @c resort. The shared staging store is reset by @ref CellStructure once it
   *  commits the staged rows (@c ensure_particle_store_synchronized), not by
   *  the decompositions -- so the type-based moves' staged rows, the two child
   *  resorts' staged rows, and the deferred commit=false hot path all keep
   * valid references until that single commit. The hybrid keeps its own copy of
   * the handle for its type-based moves. */
  void set_migration_staging(MigrationStaging staging) {
    m_regular_decomposition.set_migration_staging(staging);
    m_n_square.set_migration_staging(staging);
    m_migration_staging = std::move(staging);
  }

  auto get_cell_grid() const { return m_regular_decomposition.cell_grid; }

  auto get_cell_size() const { return m_regular_decomposition.cell_size; }

  auto get_n_square_types() const { return m_n_square_types; }

  void resort(bool global, std::vector<ParticleChange> &diff) override;

  auto get_cutoff_regular() const { return m_cutoff_regular; }

  GhostComm::HaloPlan const *halo_plan() const override { return &m_halo_plan; }

  std::span<Cell *const> local_cells() const override { return m_local_cells; }
  std::span<Cell *const> ghost_cells() const override { return m_ghost_cells; }

  Cell *particle_to_cell(Particle const &p) override {
    if (is_n_square_type(p.type())) {
      return m_n_square.particle_to_cell(p);
    }
    return m_regular_decomposition.particle_to_cell(p);
  }

  Cell const *particle_to_cell(Particle const &p) const override {
    if (is_n_square_type(p.type())) {
      return m_n_square.particle_to_cell(p);
    }
    return m_regular_decomposition.particle_to_cell(p);
  }

  Utils::Vector3d max_cutoff() const override {
    return m_n_square.max_cutoff();
  }

  Utils::Vector3d max_range() const override { return m_n_square.max_range(); }

  std::optional<BoxGeometry> minimum_image_distance() const override {
    return m_box;
  }

  BoxGeometry const &box() const override { return m_box; }

  /** @brief Count particles in child regular decompositions. */
  std::size_t count_particles_in_regular() const {
    return count_particles(m_regular_decomposition.get_local_cells());
  }

  /** @brief Count particles in child N-square decompositions. */
  std::size_t count_particles_in_n_square() const {
    return count_particles(m_n_square.get_local_cells());
  }

private:
  /**
   * @brief Build the plan-based halo plan combining the regular child's
   *        neighbors/local with the n-square child's collective section.
   */
  GhostComm::HaloPlan make_halo_plan();

  std::size_t count_particles(std::vector<Cell *> const &local_cells) const;
};
