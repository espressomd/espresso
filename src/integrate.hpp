/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#ifndef INTEGRATE_H
#define INTEGRATE_H

/** \file integrate.hpp    Molecular dynamics integrator.
 *
 *  For more information see \ref integrate.cpp "integrate.c".
*/

#define INTEG_METHOD_NPT_ISO   0
#define INTEG_METHOD_NVT       1

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Switch determining which Integrator to use. */
extern int integ_switch;

/** incremented if a Verlet update is done, aka particle resorting. */
extern int n_verlet_updates;

/** Time step for the integration. */
extern double time_step;
extern double time_step_half;
extern double time_step_squared;
extern double time_step_squared_half;

/** Old time step needed for rescaling of forces. */
extern double old_time_step;
/** Actual simulation time (only on MASTER NODE). */
extern double sim_time;
/** Maximal interaction cutoff. */
extern double max_cut;
/** Verlet list skin. */
extern double skin;

/** If non-zero, the particle data will be resorted before the next integration. */
extern int    resort_particles;
/** If non-zero, the forces will be recalculated before the next integration. */
extern int    recalc_forces;
/** Average number of integration steps the verlet list has been re
    used. */
extern double verlet_reuse;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** check sanity of integrator params */
void integrator_sanity_checks();

/** check sanity of npt params */
void integrator_npt_sanity_checks();

/** Initialize the used thermodynamic Ensemble (NVT or NPT) */
void integrate_ensemble_init();

/** integrate with velocity verlet integrator.
    \param n_steps number of steps to integrate.
    \param reuse_forces if nonzero, blindly trust
    the forces still stored with the particles for the first time step.
 */
void integrate_vv(int n_steps, int reuse_forces);

/** function that rescales all velocities on one node according to a
    new time step. */
void rescale_velocities(double scale); 

/*@}*/

#endif
