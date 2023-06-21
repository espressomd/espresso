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
#ifndef CORE_ENERGY_HPP
#define CORE_ENERGY_HPP
/** \file
 *  Energy calculation.
 *
 *  Implementation in energy.cpp.
 */

#include "Observable_stat.hpp"

#include <memory>

/** Parallel energy calculation. */
std::shared_ptr<Observable_stat> calculate_energy();

/** Calculate the total energy of the system. */
double mpi_calculate_potential_energy();

/** Helper function for @ref Observables::Energy. */
double mpi_observable_compute_energy();

/**
 * @brief Compute short-range energy of a particle.
 *
 * Iterates through particles inside cell and neighboring cells and computes
 * energy contribution for a specific particle.
 *
 * @param pid    Particle id
 * @return Non-bonded energy of the particle.
 */
double particle_short_range_energy_contribution(int pid);

#endif
#ifdef DIPOLE_FIELD_TRACKING
/** Calculate dipole fields. */
void calc_long_range_fields();
void mpi_calc_long_range_fields();
#endif
