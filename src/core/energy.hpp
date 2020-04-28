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
/** \file
 *  Energy calculation.
 *
 *  Implementation in energy.cpp.
 */

#ifndef _ENERGY_H
#define _ENERGY_H

/* include the energy files */
#include "ParticleRange.hpp"
#include "actor/ActorList.hpp"
#include "statistics.hpp"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat energy, total_energy;

extern ActorList energyActors;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** allocate energy arrays and initialize with zero */
void init_energies(Observable_stat *stat);

/** on the master node: calc energies only if necessary */
void master_energy_calc();

/** parallel energy calculation.
    @param result non-zero only on master node; will contain the cumulative over
   all nodes. */
void energy_calc(double *result, double time);

/** Calculate long range energies (P3M, ...). */
void calc_long_range_energies(const ParticleRange &particles);

/** Calculate the total energy */
double calculate_current_potential_energy_of_system();

/*@}*/

#endif
