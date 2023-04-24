/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

/** \name Integrator switches */
/**@{*/
#define INTEG_METHOD_NPT_ISO 0
#define INTEG_METHOD_NVT 1
#define INTEG_METHOD_STEEPEST_DESCENT 2
#define INTEG_METHOD_BD 3
#define INTEG_METHOD_SD 7
/**@}*/

/** \name Integrator error codes */
/**@{*/
#define INTEG_ERROR_RUNTIME -1
#define INTEG_ERROR_SIGINT -2
/**@}*/

/** \name Integrator flags */
/**@{*/
/// recalculate forces unconditionally (mostly used for timing)
#define INTEG_REUSE_FORCES_NEVER -1
/// recalculate forces if @ref recalc_forces is set
#define INTEG_REUSE_FORCES_CONDITIONALLY 0
/// do not recalculate forces (mostly when reading checkpoints with forces)
#define INTEG_REUSE_FORCES_ALWAYS 1
/**@}*/

/** Switch determining which integrator to use. */
extern int integ_switch;

/** Verlet list skin. */
extern double skin;

/** If true, the forces will be recalculated before the next integration. */
extern bool recalc_forces;

double interaction_range();

/** Check integrator parameters and incompatibilities between the integrator
 *  and the currently active thermostat(s).
 */
void integrator_sanity_checks();

/** Integrate equations of motion
 *  @param n_steps       Number of integration steps, can be zero
 *  @param reuse_forces  Decide when to re-calculate forces
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
 *  @return number of steps that have been integrated, or a negative error code
 */
int integrate(int n_steps, int reuse_forces);

int integrate_with_signal_handler(int n_steps, int reuse_forces,
                                  bool update_accumulators);

/** Get @c verlet_reuse */
double get_verlet_reuse();

/** Get time step */
double get_time_step();

/** Get simulation time */
double get_sim_time();

/** Increase simulation time (only on head node) */
void increment_sim_time(double amount);

/** Set new @ref time_step. */
void set_time_step(double value);

/** @brief Set new skin. */
void set_skin(double value);

/** @brief Set the simulation time. */
void set_time(double value);

void set_integ_switch(int value);

#endif
