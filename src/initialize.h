/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef INITIALZE_H
#define INITIALZE_H
/** \file initialize.h
    This file contains the hook procedures. These are the ones with names on_* and
    are called whenever something is changed in Espresso which might influence other parts. For example, the P3M code
    has to be reinitialized whenever the box size changes. The hooking mechanism allows to keep track of such changes.

    For this mechanism to work, two things have to be fulfilled. If some part of the code changes some property, it has
    to call the corresponding hook, i. e. on_particle_change if a particle property has been changed or
    on_short_range_ia_change, if a short ranged interaction has been changed.
    In turn procedures that depend on particle properties or the box size, should react to such changes in the
    corresponding hook procedure.
 */

#include <tcl.h>

/** \name Hook procedures
    These procedures are called if several significant changes to
    the system happen which may make a reinitialization of subsystems
    necessary. Note that all these functions are called on ALL nodes.
    if you need to do something only on the master node, check
    \ref this_node == 0. The use of the asynchronous mpi_* functions
    (e. g. mpi_bcast_parameter) on the master node is possible.
*/
/*@{*/

/** called once at the very beginning of the program start. */
int on_program_start(Tcl_Interp *interp);

/** called every time the simulation is continued/started, i. e.
    when switching from Tcl to the simulation core. */
void on_integration_start();

/** called before calculating observables, i.e. energy, pressure or the 	     integrator (forces). Initialize any methods here which are not
    initialized immediately (P3M, Maggs, etc.). */
void on_observable_calc();

/** called every time a particle property is changed via Tcl. */
void on_particle_change();

/** called every time the particles are resorted from node to node. */
void on_resort_particles();

/** called every time the coulomb parameters are changed. */
void on_coulomb_change();

/** called every time short ranged interaction parameters are changed. */
void on_short_range_ia_change();

/** called every time a constraint is changed. */
void on_constraint_change();

/** called whenever the cutoff has potentially changed. */
void on_max_cut_change();

/** called every time the cell structure is changed. */
void on_cell_structure_change();

/** called every time the NpT-integrator communicated the updated box-length. */
void on_NpT_boxl_change();

/** called every time other parameters (timestep,...) are changed. Note that
    this does not happen automatically. The callback procedure of the changed
    variable is responsible for that by calling \ref mpi_bcast_event (2).
    @param parameter is the FIELD_* identifier of the field changed.
*/
void on_parameter_change(int parameter);

/** called every time the number of particle types has changed (increased) */
void on_n_particle_types_change();

/** call this if you want to change ghost flags, e.g. wether ghosts
    have velocities or not.  This is a opt-in process, i. e. all
    features are turned off and have to be reactivated if necessary
    inside this procedure.  */
void on_ghost_flags_change();

void on_lb_params_change(int field);
void on_lbboundary_change();

/** called every time the walls for the lb fluid are changed */
void on_lbboundary_change();

/*@}*/

#endif
