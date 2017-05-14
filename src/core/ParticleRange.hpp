#ifndef CORE_PARTICLE_RANGE_HPP
#define CORE_PARTICLE_RANGE_HPP

#include "Cell.hpp"
#include "ParticleIterator.hpp"
#include "utils/Range.hpp"

using CellParticleIterator = ParticleIterator<Cell **, Particle>;
using ParticleRange = Utils::Range<CellParticleIterator>;

#endif
