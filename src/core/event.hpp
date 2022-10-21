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
#ifndef CORE_EVENT_HPP
#define CORE_EVENT_HPP
/** \file
 *  This file contains the hook procedures. These are the ones with names
 *  on_* and are called whenever something is changed in ESPResSo which
 *  might influence other parts. For example, the P3M code has to be
 *  reinitialized whenever the box size changes. The hooking mechanism
 *  allows to keep track of such changes.
 *
 *  For this mechanism to work, two things have to be fulfilled. If
 *  some part of the code changes some property, it has to call the
 *  corresponding hook, i.e. on_particle_change() if a particle
 *  property has been changed or on_short_range_ia_change(), if a short
 *  ranged interaction has been changed. In turn procedures that
 *  depend on particle properties or the box size, should react to
 *  such changes in the corresponding hook procedure.
 *
 *  Implementation in event.cpp.
 */

/** \name Hook procedures
 *  These procedures are called if several significant changes to
 *  the system happen which may make a reinitialization of subsystems
 *  necessary.
 */
/**@{*/

/** called once at the very beginning of the program start. */
void on_program_start();

/** called every time the simulation is continued/started, i.e.
 *  when switching from the script interface to the simulation core.
 */
void on_integration_start(double time_step);

/** called before calculating observables, i.e. energy, pressure or
 *  the integrator (forces). Initialize any methods here which are not
 *  initialized immediately (P3M etc.).
 */
void on_observable_calc();

/** called every time a particle property is changed via the script interface.
 */
void on_particle_change();

/** called every time the charge of a particle has changed. */
void on_particle_charge_change();

/** called every time the Coulomb parameters are changed.

all Coulomb methods have a short range part, aka near field
     correction. Even in case of switching off, we should call this,
     since the required cutoff might have reduced.
 */
void on_coulomb_change();
/** called every time the dipolar parameters are changed. */
void on_dipoles_change();

/** called every time short ranged interaction parameters are changed. */
void on_short_range_ia_change();

/** called every time a non-bonded interaction parameters are changed. */
void on_non_bonded_ia_change();

/** called every time a constraint is changed. */
void on_constraint_change();

/**
 * @brief Called when the box length has changed. This routine is relatively
 * fast, and changing the box length every time step as for example necessary
 * for NpT is more or less ok.
 *
 * @param skip_method_adaption skip the long-range methods adaptions
 */
void on_boxl_change(bool skip_method_adaption = false);

/** called every time a major change to the cell structure has happened,
 *  like the skin or grid have changed. This one is potentially slow.
 */
void on_cell_structure_change();

/** called every time the temperature changes. This one is potentially slow. */
void on_temperature_change();

/** @brief Called when the periodicity changes. Internally calls @ref
 * on_skin_change.
 */
void on_periodicity_change();

/** @brief Called when the skin is changed.
 */
void on_skin_change();

/** @brief Called when parameters of thermostats are changed.
 */
void on_thermostat_param_change();

/** @brief Called when the timestep changed.
 *  Internally calls @ref on_thermostat_param_change.
 */
void on_timestep_change();

/** @brief Called when the force cap changed.
 */
void on_forcecap_change();

/** @brief Called when the node_grid changed.
 */
void on_node_grid_change();

unsigned global_ghost_flags();

/** @brief Called when the LB boundary conditions are changed
 *  (geometry, slip velocity, or both).
 */
void on_lb_boundary_conditions_change();

/** @brief Update particles with properties depending on other particles,
 *  namely virtual sites and ICC charges.
 */
void update_dependent_particles();

/**@}*/

#endif
