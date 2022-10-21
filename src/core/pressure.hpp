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
/** \file
 *  Pressure calculation. Really similar to energy.hpp.
 */

#ifndef CORE_PRESSURE_HPP
#define CORE_PRESSURE_HPP

#include "Observable_stat.hpp"

#include <utils/Vector.hpp>

#include <memory>

/** Parallel pressure calculation from a virial expansion. */
std::shared_ptr<Observable_stat> calculate_pressure();

/** Helper function for @ref Observables::PressureTensor. */
Utils::Vector9d mpi_observable_compute_pressure_tensor();

#endif
