/*
  Copyright (C) 2010-2018 The ESPResSo project
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

/** \file
 *  Molecular dynamics integrator.
 *
 *  For more information see \ref integrate.cpp "integrate.cpp".
 */

#define INTEG_METHOD_NPT_ISO 0
#define INTEG_METHOD_NVT 1
#define INTEG_METHOD_STEEPEST_DESCENT 2

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
/** True iff the user has changed the skin setting. */
extern bool skin_set;

/** If true, the forces will be recalculated before the next integration. */
extern bool recalc_forces;
/** Average number of integration steps the Verlet list has been re
    used. */
extern double verlet_reuse;

/** Communicate signal handling to the Python interpreter */
extern bool set_py_interrupt;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** check sanity of integrator params */
void integrator_sanity_checks();

/** integrate equations of motion
    \param n_steps number of steps to integrate.
    \param reuse_forces if nonzero, blindly trust
    the forces still stored with the particles for the first time step.

    @details This function calls two hooks for propagation kernels such as
    velocity verlet, velocity verlet + npt box changes, and steepest_descent.
    One hook is called before and one after the force calculation.
    It is up to the propagation kernels to increment the simulation time.

    This function propagates the system according to the choice of integrator
    stored in @ref integ_switch. The general structure is:
    - if reuse_forces is zero, recalculate the forces based on the current
      state of the system
    - Loop over the number of simulation steps:
      -# initialization (e.g., RATTLE)
      -# First hook for propagation kernels
      -# Update dependent particles and properties (RATTLE, virtual sites)
      -# Calculate forces for the current state of the system. This includes
   forces added by the Langevin thermostat and the Lattice-Boltzmann-particle
   coupling
      -# Second hook for propagation kernels
      -# Update dependent properties (Virtual sites, RATTLE)
      -# Run single step algorithms (Lattice-Boltzmann propagation, collision
   detection, NPT update)
    - Final update of dependent properties and statistics/counters

    High-level documentation of the integration and thermostatting schemes
    can be found in doc/sphinx/system_setup.rst and /doc/sphinx/running.rst

 */
void integrate_vv(int n_steps, int reuse_forces);

/*@}*/

int python_integrate(int n_steps, bool recalc_forces, bool reuse_forces);

void integrate_set_nvt();
int integrate_set_npt_isotropic(double ext_pressure, double piston, int xdir,
                                int ydir, int zdir, bool cubic_box);

#endif
