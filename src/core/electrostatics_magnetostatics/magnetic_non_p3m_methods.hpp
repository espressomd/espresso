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
#ifndef MAG_NON_P3M_H
#define MAG_NON_P3M_H
/** \file
 *  All 3d non-P3M methods to deal with the magnetic dipoles
 *
 *  DAWAANR => Dipolar All With All And No Replica
 *   Handling of a system of dipoles where no replicas exist.
 *   Assumes minimum image convention for those axis in which the
 *   system is periodic.
 *
 *  MDDS => Magnetic Dipoles Direct Sum
 *   Calculate dipole-dipole interaction of a periodic system
 *   by explicitly summing the dipole-dipole interaction over several copies of
 *   the system.
 *   Uses spherical summation order.
 *
 */
#include "config.hpp"
#include <ParticleRange.hpp>

#ifdef DIPOLES
#include "Particle.hpp"

/** Calculate dipolar energy and/or force between two particles */
double calc_dipole_dipole_ia(Particle &p1, Particle &p2, bool force_flag);

/* =============================================================================
                  DAWAANR => DIPOLAR ALL WITH ALL AND NO REPLICA
   =============================================================================
*/

/** Core of the DAWAANR method: here you compute all the magnetic forces,
 *  torques and the magnetic energy for the whole system
 */
double dawaanr_calculations(bool force_flag, bool energy_flag,
                            ParticleRange const &particles);

/** Switch on DAWAANR magnetostatics.
 *  @return ES_ERROR, if not on a single CPU
 */
int dawaanr_set_params();

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS
   =============================================================================
*/

/** Sanity checks for the magnetic dipolar direct sum*/
int magnetic_dipolar_direct_sum_sanity_checks();

/** Core of the method: here you compute all the magnetic forces, torques and
 *  the energy for the whole system using direct sum
 */
double magnetic_dipolar_direct_sum_calculations(bool force_flag,
                                                bool energy_flag,
                                                ParticleRange const &particles);

/** Switch on direct sum magnetostatics.
 *  @param n_cut cut off for the explicit summation
 *  @return ES_ERROR, if not on a single CPU
 */
int mdds_set_params(int n_cut);

extern int Ncut_off_magnetic_dipolar_direct_sum;

#endif /*of ifdef DIPOLES  */
#endif /* of ifndef  MAG_NON_P3M_H */
