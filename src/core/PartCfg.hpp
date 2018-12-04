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

#include <boost/iterator/transform_iterator.hpp>

#include "ParticleCache.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "utils/Range.hpp"
#include "utils/SkipIterator.hpp"
#include "utils/serialization/Particle.hpp"

/**
 * @brief Proxy class that gets a particle range from
 *        from the global local_particles.
 */
class GetLocalParts {
  class SkipIfNullOrGhost {
  public:
    bool operator()(const LocalParticles::value_type &kv) const {
      return (kv.second == nullptr) or (kv.second->l.ghost);
    }
  };

  class SecondDeref {
  public:
    Particle &operator()(const LocalParticles::value_type &kv) const {
      return assert(kv.second), *(kv.second);
    }
  };

  using skip_it =
      Utils::SkipIterator<LocalParticles::iterator, SkipIfNullOrGhost>;
  using iterator = boost::transform_iterator<SecondDeref, skip_it>;
  using Range = Utils::Range<iterator>;

public:
  Range operator()() const {
    auto begin = skip_it(local_particles.begin(), local_particles.end());
    auto end = skip_it(local_particles.end(), local_particles.end());

    return {iterator(begin), iterator(end)};
  }
};

using PartCfg = ParticleCache<GetLocalParts, PositionUnfolder>;
#endif
