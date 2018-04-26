#ifndef CORE_PART_CFG_HPP
#define CORE_PART_CFG_HPP

#include <boost/iterator/indirect_iterator.hpp>

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
    bool operator()(Particle const *p_ptr) const {
      return (p_ptr == nullptr) or (p_ptr->l.ghost);
    }
  };

  using skip_it = Utils::SkipIterator<Particle **, SkipIfNullOrGhost>;
  using iterator = boost::indirect_iterator<skip_it>;
  using Range = Utils::Range<iterator>;

public:
  Range operator()() const {
    auto begin =
        skip_it(local_particles, local_particles + max_seen_particle + 1,
                SkipIfNullOrGhost());
    auto end =
        skip_it(local_particles + max_seen_particle + 1,
                local_particles + max_seen_particle + 1, SkipIfNullOrGhost());

    return Utils::make_range(make_indirect_iterator(begin),
                             make_indirect_iterator(end));
  }
};

using PartCfg = ParticleCache<GetLocalParts, PositionUnfolder>;
#endif
