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
 *  This contains a timing loop for the force calculation.
 *
 *  Implementation in tuning.cpp.
 */

#ifndef TUNING_H
#define TUNING_H

/** Measure the time for some force calculations.
 *  Actually performs \ref mpi_integrate (0).
 *  This times the force calculation without
 *  propagating the system. It therefore does
 *  not include e.g. Verlet list updates.
 *  @param int_steps  Number of integration steps.
 *  @return Time per integration in milliseconds.
 */
double time_force_calc(int int_steps);

/** Set the optimal @ref skin between @p min_skin and @p max_skin
 *  by bisection to tolerance @p tol.
 */
void tune_skin(double min_skin, double max_skin, double tol, int int_steps,
               bool adjust_max_skin);

#endif
