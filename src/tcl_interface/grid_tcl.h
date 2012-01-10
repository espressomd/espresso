/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#ifndef _GRID_TCL_H
#define _GRID_TCL_H
/** \file grid.h   Domain decomposition for parallel computing.
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different nodes (or processes
 *  or threads if you want). This grid is described in \ref
 *  node_grid. Each node has a number \ref this_node and a position
 *  \ref node_pos in that grid. Each node has also 6 nearest neighbors
 *  \ref node_neighbors which are necessary for the communication
 *  between the nodes (see also \ref ghosts.c and \ref p3m.c for more
 *  details about the communication.
 *
 *  For the 6 directions \anchor directions we have the following convention:
 *
 *  \image html directions.gif "Convention for the order of the directions"
 *
 *  The Figure illustrates the direction convetion used for arrays
 *  with 6 (e.g. \ref node_neighbors, \ref #boundary) and 3 entries
 *  (e.g \ref node_grid, \ref box_l , \ref my_left,...).
 *  
 *
 *  For more information on the domain decomposition, see \ref grid.c "grid.c". 
*/
#include <tcl.h>
#include <limits.h>
#include "utils.h"
#include "errorhandling.h"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** datafield callback for \ref #periodic. Determines wether a coordinate is pbc (default). */
int tclcallback_periodicity(Tcl_Interp *interp, void *_data);

/** datafield callback for \ref node_grid. */
int tclcallback_node_grid(Tcl_Interp *interp, void *data);

/** datafield callback for \ref #periodic. Determines wether a coordinate is pbc (default). */
int tclcallback_periodicity(Tcl_Interp *interp, void *_data);

/** datafield callback for \ref box_l. Sets the box dimensions. */
int tclcallback_box_l(Tcl_Interp *interp, void *_data);

/** changes the volume by resizing the box and isotropically adjusting the particles coordinates as well */
int tclcommand_change_volume(ClientData data, Tcl_Interp *interp, int argc, char **argv);


/*@}*/
#endif
