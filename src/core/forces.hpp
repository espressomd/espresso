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
#ifndef CORE_FORCES_HPP
#define CORE_FORCES_HPP
/** \file forces.hpp Force calculation. 
 *
 *  \todo Preprocessor switches for all forces (Default: everything is turned on).
 *  \todo Implement more flexible thermostat, %e.g. which thermostat to use.
 *
 *  For more information see forces.cpp .
 */

#include "actor/Actor.hpp"
#include "actor/ActorList.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"

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

/** Calculate forces.
 *
 *  A short list, what the function is doing:
 *  <ol>
 *  <li> Initialize forces with: \ref friction_thermo_langevin (ghost forces
 with zero).
 *  <li> Calculate bonded interaction forces:<br>
 *       Loop all local particles (not the ghosts).
 *       <ul>
 *       <li> FENE
 *       <li> ANGLE (cos bend potential)
 *       </ul>
 *  <li> Calculate non-bonded short range interaction forces:<br>
 *       Loop all \ref IA_Neighbor::vList "verlet lists" of all \ref #cells.
 *       <ul>
 *       <li> Lennard-Jones.
 *       <li> Buckingham.
 *       <li> Real space part: Coulomb.
 *       <li> Ramp.
 *       </ul>
 *  <li> Calculate long range interaction forces:<br>
 Uses <a href=P3M_calc_kspace_forces> P3M_calc_kspace_forces </a>
 *  </ol>
 */
void force_calc();

/** Check if forces are NAN 
 */
void check_forces();

/** Calculate long range forces (P3M, MMM2d...). */
void calc_long_range_forces();
/*@}*/

#endif
