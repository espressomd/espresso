/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef CORE_PART_CFG_HPP
#define CORE_PART_CFG_HPP

#include "Particle.hpp"

#include <vector>

/**
 * @brief Particle cache on the master.
 *
 * This class implements cached access to all particles in a
 * particle range on the master node.
 * This implementation fetches all particles to
 * the master on first access. Updates of the particle data are
 * triggered automatically on access. The data in the cache
 * is invalidated automatically on_particle_change, and then
 * updated on the next access.
 *
 */
class PartCfg {
  /** The particle data */
  std::vector<Particle> m_parts;
  /** State */
  bool m_valid;

public:
  using value_type = Particle;
  PartCfg() : m_valid(false) {}

  /**
   * @brief Iterator pointing to the particle with the lowest
   * id.
   *
   * Returns a random access iterator that traverses the
   * particles
   * in order of ascending id. If the cache is not up-to-date,
   * an update is triggered. This iterator stays valid as long
   * as the cache is valid. Since the cache could be invalidated
   * and updated elsewhere, iterators into the cache should not
   * be stored.
   */
  auto begin() {
    if (!m_valid)
      update();

    return m_parts.begin();
  }

  /**
   * @brief Iterator pointing past the particle with the highest
   * id.
   *
   * If the cache is not up-to-date,
   * an update is triggered.
   */
  auto end() {
    if (!m_valid)
      update();

    return m_parts.end();
  }

  /**
   * @brief Returns true if the cache is up-to-date.
   *
   * If false, particle access will trigger an update.
   */
  bool valid() const { return m_valid; }

  /**
   * @brief Invalidate the cache and free memory.
   */
  void invalidate() {
    /* Release memory */
    m_parts = std::vector<Particle>();
    /* Adjust state */
    m_valid = false;
  }

  /**
   * @brief Update particle information.
   *
   * This triggers a global update. All nodes
   * sort their particle by id, and send them
   * to the master.
   */
private:
  void update();

public:
  /** Number of particles in the config.
   */
  size_t size() {
    if (!m_valid)
      update();

    return m_parts.size();
  }

  /**
   * @brief size() == 0 ?
   */
  bool empty() {
    if (!m_valid)
      update();

    return m_parts.empty();
  }
};

#endif
