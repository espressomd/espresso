/*
 * Copyright (C) 2020 The ESPResSo project
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
#ifndef ESPRESSO_PARTICLE_DECOMPOSITION_HPP
#define ESPRESSO_PARTICLE_DECOMPOSITION_HPP

#include "BoxGeometry.hpp"
#include "Cell.hpp"
#include "LocalBox.hpp"
#include "ghosts.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

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
using ParticleChange = boost::variant<RemovedParticle, ModifiedList>;

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
   * @brief Communicator for updating ghosts from the real particles.
   */
  virtual GhostCommunicator const &exchange_ghosts_comm() const = 0;
  /**
   * @brief Communicator for force reduction.
   */
  virtual GhostCommunicator const &collect_ghost_force_comm() const = 0;

  /**
   * @brief Get pointer to local cells.
   *
   * Local cells are cells that contain particles
   * that are owned by this node.
   *
   * @return List of local cells.
   */
  virtual Utils::Span<Cell *> local_cells() = 0;

  /**
   * @brief Get pointer to local cells.
   *
   * Ghost cells are cells that contain particles
   * that are owned by different nodes but interact
   * with particles on this node.
   *
   * @return List of ghost cells.
   */
  virtual Utils::Span<Cell *> ghost_cells() = 0;

  /**
   * @brief Determine which cell a particle id belongs to.
   *
   * @param p Particle to find cell for.
   * @return Pointer to cell or nullptr if not local.
   */
  virtual Cell *particle_to_cell(Particle const &p) = 0;

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
  virtual boost::optional<BoxGeometry> minimum_image_distance() const = 0;

  virtual ~ParticleDecomposition() = default;
};

#endif
