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
#ifndef VIRTUAL_SITES_HPP
#define VIRTUAL_SITES_HPP

#include "config.hpp"
#include <serialization/Particle.hpp>

#ifdef VIRTUAL_SITES
#include "virtual_sites/VirtualSites.hpp"

/** @brief get active virtual sites implementation */
const std::shared_ptr<VirtualSites> &virtual_sites();

/** @brief Set active virtual sites implementation */
void set_virtual_sites(std::shared_ptr<VirtualSites> const &v);

#ifdef VIRTUAL_SITES_RELATIVE
int vs_relate_to(int part_num, int relate_to);

// Setup the virtual_sites_relative properties of a particle so that the given
// virtual particle will follow the given real particle Local version, expects
// both particles to be accessible through local_particles and only executes the
// changes on the virtual site locally
int local_vs_relate_to(Particle *p_current, const Particle *p_relate_to);

#endif
#endif
#endif
