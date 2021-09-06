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

#include "Observable_stat.hpp"
#include "actor/ActorList.hpp"

#include <memory>

extern ActorList energyActors;

/** Parallel energy calculation. */
std::shared_ptr<Observable_stat> calculate_energy();

/** Calculate the total energy of the system. */
double calculate_current_potential_energy_of_system();

/** Helper function for @ref Observables::Energy. */
double observable_compute_energy();

#endif
