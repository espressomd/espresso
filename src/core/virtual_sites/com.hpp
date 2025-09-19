/*
 * Copyright (C) 2010-2025 The ESPResSo project
 * Copyright (C) 2010,2011 Rudolf Weeber
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
#pragma once

/** \file
 *  This file contains routine to handle virtual sites at the center of mass of
 * a bunch of other particles (say, a molecule). Forces acting on this center of
 * mass are distributed back onto the constituents. The position/velocity/mass
 * of the virtual site at center of mass is calculated from the
 * positions/velocities/masses of many particles.
 *
 *  Virtual sites are like particles, but they will not be integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move non-virtual particles
 *  - update virtual sites
 */
#pragma once

#include "config/config.hpp"
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS

#include "BoxGeometry.hpp"
#include "cell_system/CellStructure.hpp"

#include <memory>
#include <unordered_map>
#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

void vs_com_update_particles(CellStructure &cell_structure,
                             BoxGeometry const &box_geo);
void vs_com_back_transfer_forces_and_torques(
    CellStructure &cell_structure);


#endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS