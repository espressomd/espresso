// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef INITIALZE_H
#define INITIALZE_H
/** \file initialize.h

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
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

/** called every time a particle property is changed via Tcl. */
void on_particle_change();

/** called every time the system topology changes. */
void on_topology_change();

/** called every time interaction parameters are changed. */
void on_ia_change();

/** called every time other parameters (timestep,...) are changed. Note that
    this does not happen automatically. The callback procedure of the changed
    variable is responsible for that by calling \ref mpi_bcast_event (2).
*/
void on_parameter_change();

/*@}*/

#endif
