/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_PART_CFG_HPP
#define CORE_PART_CFG_HPP

#include "ParticleCache.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "serialization/Particle.hpp"
#include <utils/SkipIterator.hpp>

#include <boost/iterator/indirect_iterator.hpp>
#include <boost/range/iterator_range.hpp>

/**
 * @brief Proxy class that gets a particle range from #local_particles.
 */
class GetLocalParts {
  class SkipIfNullOrGhost {
  public:
    bool operator()(Particle const *p_ptr) const {
      return (p_ptr == nullptr) or (p_ptr->l.ghost);
    }
  };

  using skip_it = Utils::SkipIterator<Particle **, SkipIfNullOrGhost>;
  using iterator = boost::indirect_iterator<skip_it>;
  using Range = boost::iterator_range<iterator>;

public:
  Range operator()() const {
    if (local_particles == nullptr) {
      auto begin = skip_it(nullptr, nullptr, SkipIfNullOrGhost());
      return {make_indirect_iterator(begin), make_indirect_iterator(begin)};
    }

    auto begin =
        skip_it(local_particles, local_particles + max_seen_particle + 1,
                SkipIfNullOrGhost());
    auto end =
        skip_it(local_particles + max_seen_particle + 1,
                local_particles + max_seen_particle + 1, SkipIfNullOrGhost());

    return {make_indirect_iterator(begin), make_indirect_iterator(end)};
  }
};

/** @brief Cache of particles */
using PartCfg = ParticleCache<GetLocalParts, PositionUnfolder>;
#endif
