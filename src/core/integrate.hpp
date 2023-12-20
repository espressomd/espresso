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
 *  Molecular dynamics integrator.
 *
 *  Implementation in \ref integrate.cpp.
 */

#include "config/config.hpp"

#include "PropagationMode.hpp"

#ifdef WALBERLA
#include <string>

#include <utils/Vector.hpp>
#endif

/** \name Integrator error codes */
/**@{*/
#define INTEG_ERROR_RUNTIME -1
#define INTEG_ERROR_SIGINT -2
/**@}*/

/** \name Integrator flags */
/**@{*/
/// recalculate forces unconditionally (mostly used for timing)
#define INTEG_REUSE_FORCES_NEVER -1
/// recalculate forces only if @ref Propagation::recalc_forces is set
#define INTEG_REUSE_FORCES_CONDITIONALLY 0
/// do not recalculate forces (mostly when reading checkpoints with forces)
#define INTEG_REUSE_FORCES_ALWAYS 1
/**@}*/

#ifdef WALBERLA
void walberla_tau_sanity_checks(std::string method, double tau,
                                double time_step);
void walberla_agrid_sanity_checks(std::string method,
                                  Utils::Vector3d const &geo_left,
                                  Utils::Vector3d const &geo_right,
                                  Utils::Vector3d const &lattice_left,
                                  Utils::Vector3d const &lattice_right,
                                  double agrid);
#endif // WALBERLA
