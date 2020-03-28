/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef CORE_PARTICLE_CACHE_HPP
#define CORE_PARTICLE_CACHE_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <utility>

#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/NoOp.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/serialization/flat_set.hpp>

#include <boost/algorithm/clamp.hpp>

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
 * To update the cache particles are sorted by id on the nodes,
 * and the sorted arrays a merged in a reduction tree, until the
 * master node receives a complete and sorted particle array.
 *
 * This class can be customized by running a unary operation on
 * the particles. This op is run on all the nodes. It can be used
 * e.g. to fold or unfold the coordinates on the fly.
 *
 * To iterate over the particles using the iterators is more
 * efficient than using operator[].
 *
 */
class ParticleCache {
  /** The particle data */
  std::vector<Particle> remote_parts;
  /** State */
  bool m_valid;

public:
  using value_type = Particle;
  ParticleCache() : m_valid(false) {}

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

    return remote_parts.begin();
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

    return remote_parts.end();
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
    remote_parts = std::vector<Particle>();
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
  void update() {
    if (m_valid)
      return;

    remote_parts.clear();

    auto const ids = get_particle_ids();
    auto const chunk_size = fetch_cache_max_size();

    for (size_t offset = 0; offset < ids.size();) {
      auto const this_size =
          boost::algorithm::clamp(chunk_size, 0, ids.size() - offset);
      auto const chunk_ids =
          Utils::make_const_span(ids.data() + offset, this_size);

      prefetch_particle_data(chunk_ids);

      for (auto id : chunk_ids) {
        remote_parts.push_back(get_particle_data(id));

        auto &p = remote_parts.back();
        p.r.p += image_shift(p.l.i, box_geo.length());
        p.l.i = {};
      }

      offset += this_size;
    }

    m_valid = true;
  }

  /** Number of particles in the config.
   */
  size_t size() {
    if (!m_valid)
      update();

    return remote_parts.size();
  }

  /**
   * @brief size() == 0 ?
   */
  bool empty() {
    if (!m_valid)
      update();

    return remote_parts.empty();
  }
};

#endif
