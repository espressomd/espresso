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
#ifndef VIRTUAL_SITES_HPP
#define VIRTUAL_SITES_HPP

#include "config.hpp"

#ifdef VIRTUAL_SITES
#include "Particle.hpp"

#include "virtual_sites/VirtualSites.hpp"

#include <memory>

/** @brief get active virtual sites implementation */
const std::shared_ptr<VirtualSites> &virtual_sites();

/** @brief Set active virtual sites implementation */
void set_virtual_sites(std::shared_ptr<VirtualSites> const &v);

#ifdef VIRTUAL_SITES_RELATIVE

/** Setup the @ref ParticleProperties::vs_relative "vs_relative" of a particle
 *  so that the given virtual particle will follow the given real particle.
 */
void vs_relate_to(int part_num, int relate_to);

/** Setup the @ref ParticleProperties::vs_relative "vs_relative" of a particle
 *  so that the given virtual particle will follow the given real particle.
 */
void local_vs_relate_to(Particle &p_current, Particle const &p_relate_to);

#endif
#endif
#endif
