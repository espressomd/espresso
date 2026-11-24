/*
 * Copyright (C) 2020-2026 The ESPResSo project
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

#include "cell_system/Cell.hpp"

#include "BoxGeometry.hpp"
#include "ghosts/HaloPlan.hpp"

#include <utils/Vector.hpp>

#include <cassert>
#include <optional>
#include <span>
#include <variant>
#include <vector>

struct RemovedParticle {
  int id;
};

struct ModifiedList {
  ParticleList &pl;
};

/**
 * @brief Change of Particle Address.
 */
using ParticleChange = std::variant<RemovedParticle, ModifiedList>;

/**
 * @brief Data kept outside of the particles that has to follow them when
 * they move to another MPI rank, e.g. their bonds.
 *
 * A particle decomposition calls @ref pack() for every particle it sends
 * to another rank, in the order of the particles in the message, and
 * @ref unpack() for every received particle, in the same order.
 */
class ParticleMigration {
public:
  /** @brief Move the data of a departing particle into a buffer. */
  virtual void pack(int pid, std::vector<int> &buffer) = 0;
  /**
   * @brief Restore the data of an arriving particle from a buffer.
   * @return the remainder of the buffer.
   */
  virtual std::span<const int> unpack(int pid, std::span<const int> buffer) = 0;
  virtual ~ParticleMigration() = default;
};

/**
 * @brief A distributed particle decomposition.
 *
 * An implementation of this class organizes particles into cells.
 * It owns the particles, and provides a way of restoring the order
 * when it is disturbed, and provides a description of the neighborhood
 * relations between the cells, by which pair interactions with particles
 * near-by can be calculated. Related to this it provides descriptions
 * of the ghost communications by which particles can be synchronized that
 * are not owned locally, but interact with local particles.
 */
class ParticleDecomposition {
public:
  /**
   * @brief Resort particles.
   *
   * After calling this function, every particle is in its home cell.
   * The output parameter is filled with the changes to the local
   * particle content, which allows e.g. to keep particles indices
   * in an efficient way.
   *
   * This is a collective call.
   *
   * @param[in] global_flag Expect particles to be displaced by more than a
   * local box size.
   * @param[out] diff Cells that have been touched.
   */
  virtual void resort(bool global_flag, std::vector<ParticleChange> &diff) = 0;

  /**
   * @brief Topology-agnostic halo-exchange plan.
   *
   * Returns a @ref GhostComm::HaloPlan describing the ghost communication for
   * this decomposition. All decompositions implement this.
   */
  virtual GhostComm::HaloPlan const *halo_plan() const = 0;

  /**
   * @brief Get pointer to local cells.
   *
   * Local cells are cells that contain particles
   * that are owned by this node.
   *
   * @return List of local cells.
   */
  virtual std::span<Cell *const> local_cells() const = 0;

  /**
   * @brief Get pointer to local cells.
   *
   * Ghost cells are cells that contain particles
   * that are owned by different nodes but interact
   * with particles on this node.
   *
   * @return List of ghost cells.
   */
  virtual std::span<Cell *const> ghost_cells() const = 0;

  /**
   * @brief Determine which cell a particle id belongs to.
   *
   * @param p Particle to find cell for.
   * @return Pointer to cell or nullptr if not local.
   */
  virtual Cell *particle_to_cell(Particle const &p) = 0;
  virtual Cell const *particle_to_cell(Particle const &p) const = 0;

  /**
   * @brief Maximum supported cutoff.
   */
  virtual Utils::Vector3d max_cutoff() const = 0;

  /**
   * @brief Range in which calculations are performed.
   */
  virtual Utils::Vector3d max_range() const = 0;

  /**
   * @brief Return the box geometry needed for distance calculation
   *        if minimum image convention should be used needed for
   *        distance calculation.
   */
  virtual std::optional<BoxGeometry> minimum_image_distance() const = 0;

  virtual BoxGeometry const &box() const = 0;

  /**
   * @brief Set the data that follows the particles when they move to
   * another MPI rank. Must be set identically on all ranks.
   */
  virtual void set_migration(ParticleMigration *migration) {
    m_migration = migration;
  }

  virtual ~ParticleDecomposition() = default;

protected:
  ParticleMigration *m_migration = nullptr;

  /** @brief Pack the migration data of the particles of a container. */
  template <class Particles>
  void pack_migration(Particles const &particles, std::vector<int> &buffer) {
    buffer.clear();
    if (m_migration) {
      for (auto const &p : particles) {
        m_migration->pack(p.id(), buffer);
      }
    }
  }

  /** @brief Unpack the migration data of the particles of a container. */
  template <class Particles>
  void unpack_migration(Particles const &particles,
                        std::vector<int> const &buffer) {
    if (m_migration) {
      auto rest = std::span<const int>(buffer);
      for (auto const &p : particles) {
        rest = m_migration->unpack(p.id(), rest);
      }
      assert(rest.empty());
    }
  }
};
