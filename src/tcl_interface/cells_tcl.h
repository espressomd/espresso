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
#ifndef _CELLS_TCL_H
#define _CELLS_TCL_H
/** \file cells.h
    This file contains everything related to the cell structure / cell
    system.
    
    The cell system (\ref Cell Structure) describes how particles are
    distributed on the cells and how particles of different cells
    (regardless if they reside on the same or different nodes)
    interact with each other. The following cell systems are implemented:
  
    <ul>
    <li> domain decomposition: The simulation box is divided spatially
    ino cells (see \ref domain_decomposition.h). This is suitable for
    short range interctions.
    <li> nsquare: The particles are distributed equally on all nodes
    regardless their spatial position (see \ref nsquare.h). This is
    suitable for long range interactions that can not be treated by a
    special method like P3M (see \ref p3m.h).
    <li> layered: in x and y directions, it uses a nsquared type of interaction calculation,
                  but in z it has a domain decomposition into layers.
    </ul>
  
    One can switch between different cell systems with the tcl command
    cellsystem implemented in \ref cells.c .
  
    Some structures are common to all cell systems: 
  
   <ul>
   <li> All cells, real cells as well as ghost cells, are stored in the array \ref cells::cells with size \ref
   n_cells. The size of this array has to be changed with \ref realloc_cells.
   <li> Their are two lists of cell pointers to acces particles and
   ghost particles on a node: \ref local_cells contains pointers to
   all cells containing the particles physically residing on that
   node. \ref ghost_cells contains pointers to all cells containing
   the ghost particles of that node. The size of these lists has to be
   changed with \ref realloc_cellplist
   <li> An example using the cell pointer lists to access particle data
   can be found in the function \ref
   print_local_particle_positions. DO NOT INVENT YOUR OWN WAY!!!
   </ul>
*/

#include "particle_data.h"


/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** implementation of the Tcl command cellsystem */
int tclcommand_cellsystem(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv);

#endif
