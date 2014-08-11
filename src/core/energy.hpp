/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
/** \file energy.hpp
    Implementation of the energy calculation.
*/
/*
#include "utils.hpp"
#include "integrate.hpp"
#include "object-in-fluid/stretching_force.hpp"
#include "object-in-fluid/stretchlin_force.hpp"
#include "object-in-fluid/area_force_local.hpp"
#include "object-in-fluid/area_force_global.hpp"
#include "object-in-fluid/bending_force.hpp"
#include "object-in-fluid/volume_force.hpp"
#include "dihedral.hpp"
#include "mdlc_correction.hpp"
*/
#ifndef _ENERGY_H
#define _ENERGY_H

/* include the energy files */
#include "statistics.hpp"
#include "actor/ActorList.hpp"

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
    @param result non-zero only on master node; will contain the cumulative over all nodes. */
void energy_calc(double *result);

/** Calculate long range energies (P3M, MMM2d...). */
void calc_long_range_energies();

/* ************************** inline ************************** */

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param ia_params the interaction parameters between the two particles
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2.
    @return the short ranged interaction energy between the two particles
*/
double calc_non_bonded_pair_energy(Particle *p1, Particle *p2,
					    IA_parameters *ia_params,
                        double d[3], double dist, double dist2);

/** Add non bonded energies and short range coulomb between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
void add_non_bonded_pair_energy(Particle *p1, Particle *p2, double d[3],
                     double dist, double dist2);

/** Calculate bonded energies for one particle.
    @param p1 particle for which to calculate energies
*/
void add_bonded_energy(Particle *p1);

/** Calculate kinetic energies for one particle.
    @param p1 particle for which to calculate energies
*/
void add_kinetic_energy(Particle *p1);

/*@}*/

#endif
