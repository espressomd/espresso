/*
 * Copyright (C) 2012-2019 The ESPResSo project
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
#ifndef _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_H
/** \file
 *  Routines to calculate the OIF global forces energy or/and and force
 *  for a particle triple (triangle from mesh). See @cite dupin07a.
 */

#include "CellStructure.hpp"
#include "oif_global_forces_params.hpp"

#include <utils/Vector.hpp>

/** Calculate the OIF global force.
 *  Called in force_calc() from within forces.cpp
 *  - calculates the global area and global volume for a cell before the forces
 *    are handled
 *  - MPI synchronization with all reduce
 *  - !!! loop over particles from domain_decomposition !!!
 */
Utils::Vector2d calc_oif_global(int molType, CellStructure &cs);

/** Distribute the OIF global forces to all particles in the mesh. */
void add_oif_global_forces(Utils::Vector2d const &area_volume, int molType,
                           CellStructure &cs);

extern int max_oif_objects;

void mpi_set_max_oif_objects(int max_oif_objects);
#endif
