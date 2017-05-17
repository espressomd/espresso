#ifndef CORE_PART_CFG_HPP
#define CORE_PART_CFG_HPP

#include "cells.hpp"
#include "ParticleCache.hpp"

/** Particles' current configuration. Before using that
    call \ref updatePartCfg or \ref sortPartCfg to allocate
    the data if necessary (which is decided by \ref updatePartCfg). */
extern ParticleCache<CellPList> partCfg;

#endif
