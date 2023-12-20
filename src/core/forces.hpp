/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#pragma once

/** \file
 *  Force calculation.
 *
 *  Implementation in forces.cpp.
 */

#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

/** Assign external forces/torques to real particles and zero to ghosts. */
void init_forces(ParticleRange const &particles, double time_step);

/** Set forces of all ghosts to zero */
void init_forces_ghosts(ParticleRange const &particles);

/** Calculate long range forces (P3M, ...). */
void calc_long_range_forces(ParticleRange const &particles);

#ifdef NPT
/** Update the NpT virial */
void npt_add_virial_force_contribution(Utils::Vector3d const &force,
                                       Utils::Vector3d const &d);
#endif
