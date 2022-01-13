/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef CORE_FORCES_HPP
#define CORE_FORCES_HPP
/** \file
 *  Force calculation.
 *
 *  \todo Preprocessor switches for all forces (Default: everything is
 *        turned on).
 *  \todo Implement more flexible thermostat, e.g. which thermostat to use.
 *
 *  Implementation in forces.cpp.
 */

#include "CellStructure.hpp"
#include "ParticleRange.hpp"
#include "actor/ActorList.hpp"

#include <utils/Vector.hpp>

extern ActorList forceActors;

/** initialize real particle forces with thermostat forces and
    ghost particle forces with zero. */
void init_forces(const ParticleRange &particles, double time_step);

/** Set forces of all ghosts to zero */
void init_forces_ghosts(const ParticleRange &particles);

/** Calculate forces.
 *
 *  A short list, what the function is doing:
 *  <ol>
 *  <li> Initialize forces
 *  <li> Calculate bonded interaction forces
 *  <li> Calculate non-bonded short range interaction forces
 *  <li> Calculate long range interaction forces
 *  </ol>
 */
void force_calc(CellStructure &cell_structure, double time_step, double kT);

/** Calculate long range forces (P3M, ...). */
void calc_long_range_forces(const ParticleRange &particles);

#ifdef NPT
/** Update the NpT virial */
void npt_add_virial_force_contribution(const Utils::Vector3d &force,
                                       const Utils::Vector3d &d);
#endif

#endif
