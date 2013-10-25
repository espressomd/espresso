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
#ifndef MAG_NON_P3M_H
#define MAG_NON_P3M_H
/** \file magnetic_non_p3m_methods.hpp   Header of all 3d non P3M methods to deal with the magnetic dipoles
 *   
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *  Handling of a system of dipoles where no replicas exist
 *  assuming minimum image convention
 *
 *  MDDS => Magnetic dipoles direct sum, compute the interactions via direct sum, 
 *
 */
#include "utils.hpp"

#ifdef DIPOLES
#include "particle_data.hpp"

// Calculates dipolar energy and/or force between two particles
double calc_dipole_dipole_ia(Particle* p1, Particle *p2, int force_flag);

/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/

/** Core of the DAWAANR method: here you compute all the magnetic forces, torques and the magnetic energy for the whole system*/
double dawaanr_calculations(int force_flag, int energy_flag) ;

/** switch on DAWAANR magnetostatics.
    @return ES_ERROR, if not on a single CPU
 */
int dawaanr_set_params();

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS               
   =============================================================================
*/

/* Sanity checks for the magnetic dipolar direct sum*/
int magnetic_dipolar_direct_sum_sanity_checks();

/* Core of the method: here you compute all the magnetic forces,torques and the energy for the whole system using direct sum*/
double  magnetic_dipolar_direct_sum_calculations(int force_flag, int energy_flag);

/** switch on direct sum magnetostatics.
    @param n_cut cut off for the explicit summation
    @return ES_ERROR, if not on a single CPU
 */
int mdds_set_params(int n_cut);

#endif /*of ifdef DIPOLES  */
#endif /* of ifndef  MAG_NON_P3M_H */
