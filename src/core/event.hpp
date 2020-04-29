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
 *  necessary. Note that all these functions are called on ALL nodes.
 *  If you need to do something only on the master node, check
 *  \ref this_node == 0. The use of the asynchronous mpi_* functions
 *  (e.g. mpi_bcast_parameter) on the master node is possible.
 */
/*@{*/

/** called once at the very beginning of the program start. */
void on_program_start();

/** called every time the simulation is continued/started, i.e.
 *  when switching from the script interface to the simulation core.
 */
void on_integration_start();

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

/** called every time the particles are resorted from node to node. */
void on_resort_particles();

/** called every time the Coulomb parameters are changed. */
void on_coulomb_change();

/** called every time short ranged interaction parameters are changed. */
void on_short_range_ia_change();

/** called every time a constraint is changed. */
void on_constraint_change();

/** called every time the box length has changed. This routine
 *  is relatively fast, and changing the box length every time step
 *  as for example necessary for NpT is more or less ok.
 */
void on_boxl_change();

/** called every time a major change to the cell structure has happened,
 *  like the skin or grid have changed. This one is potentially slow.
 */
void on_cell_structure_change();

/** called every time the temperature changes. This one is potentially slow. */
void on_temperature_change();

/** called every time other parameters (timestep,...) are changed. Note that
 *  this does not happen automatically. The callback procedure of the changed
 *  variable is responsible for that.
 *  @param parameter is the @ref Fields identifier of the field changed.
 */
void on_parameter_change(int parameter);

unsigned global_ghost_flags();

/** called every time the walls for the lb fluid are changed */
void on_lbboundary_change();

/** @brief Update particles with properties depending on other particles
 *   e.g., virtual sites, ICC charges
 */
void update_dependent_particles();

/*@}*/

#endif
