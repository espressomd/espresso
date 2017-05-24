#ifndef CORE_PART_CFG_HPP
#define CORE_PART_CFG_HPP

#include "cells.hpp"
#include "grid.hpp"
#include "ParticleCache.hpp"
#include "utils/serialization/Particle.hpp"

/** Particles' current configuration. Particle coordinates
 are unfolded. */
extern ParticleCache<CellPList, PositionUnfolder> partCfg;

#endif
