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
#ifndef ESPRESSO_SRC_CORE_DPD_HPP
#define ESPRESSO_SRC_CORE_DPD_HPP
/** \file
 *  Routines to use DPD as thermostat or pair force @cite soddemann03a
 *
 *  Implementation in @ref dpd.cpp.
 */

#include "config/config.hpp"

#ifdef DPD

#include "Particle.hpp"

#include <utils/Vector.hpp>

struct IA_parameters;

void dpd_init(double kT, double time_step);

Utils::Vector3d dpd_pair_force(Particle const &p1, Particle const &p2,
                               IA_parameters const &ia_params,
                               Utils::Vector3d const &d, double dist,
                               double dist2);
Utils::Vector9d dpd_stress();

#endif // DPD
#endif
