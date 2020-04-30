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
#ifndef INTEGRATE_H
#define INTEGRATE_H

/** \file
 *  Molecular dynamics integrator.
 *
 *  Implementation in \ref integrate.cpp.
 */

#define INTEG_METHOD_NPT_ISO 0
#define INTEG_METHOD_NVT 1
#define INTEG_METHOD_STEEPEST_DESCENT 2
#define INTEG_METHOD_BD 3

/** Switch determining which integrator to use. */
extern int integ_switch;

/** incremented if a Verlet update is done, aka particle resorting. */
extern int n_verlet_updates;

/** Time step for the integration. */
extern double time_step;
extern double time_step_half;
extern double time_step_squared;
extern double time_step_squared_half;

/** Actual simulation time (only on MASTER NODE). */
extern double sim_time;
/** Verlet list skin. */
extern double skin;
/** True iff the user has changed the skin setting. */
extern bool skin_set;

/** If true, the forces will be recalculated before the next integration. */
extern bool recalc_forces;
/** Average number of integration steps the Verlet list has been re-using. */
extern double verlet_reuse;

/** Communicate signal handling to the Python interpreter */
extern bool set_py_interrupt;

double interaction_range();

/** check sanity of integrator params */
void integrator_sanity_checks();

/** Integrate equations of motion
 *  @param n_steps       Number of integration steps, can be zero
 *  @param reuse_forces  Decide when to re-calculate forces:
 *                       - -1: recalculate forces unconditionally
 *                         (mostly used for timing)
 *                       - 0: recalculate forces if @ref recalc_forces is set,
 *                         meaning it is probably necessary
 *                       - 1: do not recalculate forces (mostly when reading
 *                         checkpoints with forces)
 *
 *  @details This function calls two hooks for propagation kernels such as
 *  velocity verlet, velocity verlet + npt box changes, and steepest_descent.
 *  One hook is called before and one after the force calculation.
 *  It is up to the propagation kernels to increment the simulation time.
 *
 *  This function propagates the system according to the choice of integrator
 *  stored in @ref integ_switch. The general structure is:
 *  - if reuse_forces is zero, recalculate the forces based on the current
 *    state of the system
 *  - Loop over the number of simulation steps:
 *    -# initialization (e.g., RATTLE)
 *    -# First hook for propagation kernels
 *    -# Update dependent particles and properties (RATTLE, virtual sites)
 *    -# Calculate forces for the current state of the system. This includes
 *       forces added by the Langevin thermostat and the
 *       Lattice-Boltzmann-particle coupling
 *    -# Second hook for propagation kernels
 *    -# Update dependent properties (Virtual sites, RATTLE)
 *    -# Run single step algorithms (Lattice-Boltzmann propagation, collision
 *       detection, NpT update)
 *  - Final update of dependent properties and statistics/counters
 *
 *  High-level documentation of the integration and thermostatting schemes
 *  can be found in doc/sphinx/system_setup.rst and /doc/sphinx/running.rst
 *
 *  @return number of steps that have been integrated
 */
int integrate(int n_steps, int reuse_forces);

/** @brief Run the integration loop. Can be interrupted with Ctrl+C.
 *
 *  @param n_steps        Number of integration steps, can be zero
 *  @param recalc_forces  Whether to recalculate forces
 *  @param reuse_forces   Whether to re-use forces
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int python_integrate(int n_steps, bool recalc_forces, bool reuse_forces);

/** @brief Set the steepest descent integrator for energy minimization.
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int integrate_set_steepest_descent(double f_max, double gamma, int max_steps,
                                   double max_displacement);

/** @brief Set the velocity Verlet integrator for the NVT ensemble. */
void integrate_set_nvt();

/** @brief Set the Brownian Dynamics integrator. */
void integrate_set_bd();

/** @brief Set the velocity Verlet integrator modified for the NpT ensemble
 *  with isotropic rescaling.
 *
 *  @param ext_pressure  Reference pressure
 *  @param piston        Piston mass
 *  @param xdir_rescale  Enable box rescaling in the *x*-direction
 *  @param ydir_rescale  Enable box rescaling in the *y*-direction
 *  @param zdir_rescale  Enable box rescaling in the *z*-direction
 *  @param cubic_box     Determines if the volume fluctuations should also
 *                       apply to dimensions which are switched off by the
 *                       above flags and which do not contribute to the
 *                       pressure (3D) or tension (2D, 1D)
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int integrate_set_npt_isotropic(double ext_pressure, double piston,
                                bool xdir_rescale, bool ydir_rescale,
                                bool zdir_rescale, bool cubic_box);

#endif
