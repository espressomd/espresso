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
 *  for a particle triple (triangle from mesh). (Dupin2007)
 */

#include <ParticleRange.hpp>

/** Set parameters for the OIF global forces potential. */
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g,
                                 double V0, double kv);

/** Calculate the OIF global force.
 *  Called in force_calc() from within forces.cpp
 *  - calculates the global area and global volume for a cell before the forces
 *    are handled
 *  - sums up parts for area with mpi_reduce from local triangles
 *  - synchronization with allreduce
 *  - !!! loop over particles from domain_decomposition !!!
 */
void calc_oif_global(double *area_volume, int molType,
                     const ParticleRange &particles);

/** Distribute the OIF global forces to all particles in the mesh. */
void add_oif_global_forces(double const *area_volume, int molType,
                           const ParticleRange &particles);

extern int max_oif_objects;
#endif
