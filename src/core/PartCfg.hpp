/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "Particle.hpp"

#include <vector>

/**
 * @brief Particle cache on the head node.
 *
 * This class implements cached access to all particles in a
 * particle range on the head node.
 * This implementation fetches all particles to the head node on creation.
 */
class PartCfg {
  /** The particle data */
  std::vector<Particle> m_parts;
  BoxGeometry const &m_box_geo;

public:
  using value_type = Particle;
  explicit PartCfg(BoxGeometry const &box_geo) : m_parts{}, m_box_geo{box_geo} {
    update();
  }

  /** @brief Iterator pointing to the particle with the lowest id. */
  auto begin() { return m_parts.begin(); }

  /** @brief Iterator pointing past the particle with the highest id. */
  auto end() { return m_parts.end(); }

  /** @brief Number of particles in the config. */
  auto size() { return m_parts.size(); }

  /** @brief Is the config empty? */
  auto empty() { return m_parts.empty(); }

private:
  /**
   * @brief Update particle information.
   *
   * This triggers a global update. All nodes
   * sort their particle by id, and send them
   * to the head node.
   */
  void update();
};
