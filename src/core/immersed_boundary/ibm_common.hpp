/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#ifndef IMMERSED_BOUNDARY_IBM_COMMON_HPP
#define IMMERSED_BOUNDARY_IBM_COMMON_HPP

#include <utils/Vector.hpp>

/**
 * @brief Returns the position of a given particle.
 *
 * @param pid %Particle id.
 * @return position of the particle.
 */
Utils::Vector3d get_ibm_particle_position(int pid);

#endif
