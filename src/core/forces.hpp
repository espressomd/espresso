/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _FORCES_HPP
#define _FORCES_HPP
/** \file forces.hpp Force calculation. 
 *
 *  \todo Preprocessor switches for all forces (Default: everything is turned on).
 *  \todo Implement more flexible thermostat, %e.g. which thermostat to use.
 *
 *  For more information see forces.cpp .
 */

#include "iccp3m.hpp"
#include "external_potential.hpp"
#include "actor/Actor.hpp"
#include "actor/ActorList.hpp"
extern ActorList forceActors;

/** \name Exported Functions */
/************************************************************/
/*@{*/

/******************* forces.cpp *******************/

/** initialize real particle forces with thermostat forces and
    ghost particle forces with zero. */
void init_forces();

/** Set forces of all ghosts to zero
 */
void init_forces_ghosts();

/** Check if forces are NAN 
 */
void check_forces();

/** Calculate long range forces (P3M, MMM2d...). */
void calc_long_range_forces();

void 
calc_non_bonded_pair_force_from_partcfg(Particle *p1, Particle *p2, 
                                        IA_parameters *ia_params,
                                        double d[3], double dist, double dist2,
                                        double force[3],
                                        double torque1[3] = NULL, 
                                        double torque2[3] = NULL);

void
calc_non_bonded_pair_force_from_partcfg_simple(Particle *p1, Particle *p2,
                                               double d[3], double dist,
                                               double dist2, double force[3]);
/*@}*/

#endif
