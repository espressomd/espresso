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
#include "serialization/Particle.hpp"

/**
 * @brief Proxy class that gets a particle range from
 *        from the global local_particles.
 */
class GetLocalParts {
public:
  ParticleRange operator()() const { return local_cells.particles(); }
};

using PartCfg = ParticleCache<GetLocalParts, PositionUnfolder>;
#endif
