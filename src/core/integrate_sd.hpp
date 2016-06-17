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
#ifndef __INTEGRATE_SD_H
#define __INTEGRATE_SD_H

/** Average number of integration steps the verlet list has been re
    used. */
extern double verlet_reuse;

extern double sd_viscosity;

extern double sd_radius;

extern int sd_seed[2];
extern int sd_random_state[2];

extern double sd_random_precision;

#ifdef SD

#ifdef CUDA
#  define HAVE_CUBLAS
#  include <cublas_v2.h>
//#include <magma.h>
#  ifdef __CUDA_ARCH__
#    if __CUDA_ARCH__ < 200
#      error "Stokesian Dynamics needs at least CUDA ComputeCapability 2.0"
#    endif //#if __CUDA_ARCH__ < 200
#  endif //#ifdef  __CUDA_ARCH__
#else
#  ifdef SD
#    error "CUDA is not given!"
#    error "StokesDynamics requires CUDA"
#  endif
#endif

#endif

#define DIM (3)

/** \file integrate_sd.hpp    Stokes dynamics integrator.
 *
 *  For more information see \ref integrate_sd.cpp "integrate_sd.cpp".
*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

#ifdef SD_USE_FLOAT
typedef float real;
#else
typedef double real;
#endif

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

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** check sanity of integrator params */
void integrator_sanity_checks_sd();

/** integrate with euler integrator.
    \param n_steps number of steps to integrate.
 */
void integrate_sd(int n_steps);

/** set particle apart so they dont overlap.
 */
int sd_set_particles_apart();

int sd_test(int size, int type);

/*@}*/
//#endif /* SD */
#endif
