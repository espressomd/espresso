/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file statistics_fluid.hpp
 *
 * Fluid related analysis functions.
 * Header file for \ref statistics_fluid.cpp.
 *
 */

#ifndef STATISTICS_FLUID_H
#define STATISTICS_FLUID_H

#include "utils.hpp"

#ifdef LB

/** Caclulate mass of the LB fluid.
 * \param result Fluid mass
 */
void lb_calc_fluid_mass(double *result);

/** Calculate momentum of the LB fluid.
 * \param result Fluid momentum
 */
void lb_calc_fluid_momentum(double *result);

/** Calculate temperature of the LB fluid.
 * \param result Fluid temperature
 */
void lb_calc_fluid_temp(double *result);

void lb_collect_boundary_forces(double *result);

void lb_calc_densprof(double *result, int *params);
void lb_calc_velprof(double *result, int *params);

#endif /* LB */

#endif /* STATISTICS_FLUID_H */
