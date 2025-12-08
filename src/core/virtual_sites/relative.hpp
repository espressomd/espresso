/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

/** Get real particle tracked by a virtual site.
 *  @param cell_structure Cell structure.
 *  @param p Virtual site.
 *  @return Pointer to real particle, or nullptr if lookup fails.
 */
Particle *get_reference_particle(CellStructure &cell_structure,
                                 Particle const &p);

void vs_relative_update_particles(CellStructure &cell_structure,
                                  BoxGeometry const &box_geo);
void vs_relative_back_transfer_forces_and_torques(
    CellStructure &cell_structure);
Utils::Matrix<double, 3, 3>
vs_relative_pressure_tensor(CellStructure const &cell_structure);

#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
