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
#ifndef RATTLE_H
#define RATTLE_H

/** \file
 *  RATTLE algorithm (@cite andersen83a).
 *
 *  For more information see \ref rattle.cpp.
 */

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/rigid_bond.hpp"
#include "cells.hpp"
#include "config.hpp"

#ifdef BOND_CONSTRAINT

/** Transfer the current particle positions from @ref ParticlePosition::p
 *  "Particle::r::p" to @ref ParticlePosition::p_old "Particle::r::p_old"
 */
void save_old_pos(const ParticleRange &particles,
                  const ParticleRange &ghost_particles);

/** Propagate velocity and position while using SHAKE algorithm for bond
 *  constraint.
 */
void correct_pos_shake(CellStructure &cs);

/** Correction of current velocities using RATTLE algorithm. */
void correct_vel_shake(CellStructure &cs);

#endif
#endif
