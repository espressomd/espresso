/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef VIRTUAL_SITES_RELATIVE

#include "BoxGeometry.hpp"
#include "cell_system/CellStructure.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

void vs_relative_update_particles(CellStructure &cell_structure,
                                  BoxGeometry const &box_geo);
void vs_relative_back_transfer_forces_and_torques(
    CellStructure &cell_structure);
Utils::Matrix<double, 3, 3>
vs_relative_pressure_tensor(CellStructure const &cell_structure);

#endif // VIRTUAL_SITES_RELATIVE
