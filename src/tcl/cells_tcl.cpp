/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file cells.cpp
 *
 *  This file contains functions for the cell system.
 *
 *  For more information on cells, see cells.hpp
 *   */
#include "parser.hpp"
#include "domain_decomposition.hpp"
#include "layered.hpp"
#include "ghosts.hpp"
#include "verlet.hpp"

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

int tclcommand_sort_particles(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  mpi_bcast_event(SORT_PARTICLES);
  return TCL_OK;
}

int tclcommand_cellsystem(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv)
{
  int err = 0;

  if (argc <= 1) {
    Tcl_AppendResult(interp, "usage: cellsystem <system> <params>", (char *)NULL);
    return TCL_ERROR;
  }

  if (ARG1_IS_S("domain_decomposition")) {
    if (argc > 2) {
      if (ARG_IS_S(2,"-verlet_list"))
	dd.use_vList = 1;
      else if(ARG_IS_S(2,"-no_verlet_list")) 
	dd.use_vList = 0;
      else{
	Tcl_AppendResult(interp, "wrong flag to",argv[0],
			 " : should be \" -verlet_list or -no_verlet_list \"",
			 (char *) NULL);
	return (TCL_ERROR);
      }
    }
    /** by default use verlet list */
    else dd.use_vList = 1;
    mpi_bcast_cell_structure(CELL_STRUCTURE_DOMDEC);
  }
  else if (ARG1_IS_S("nsquare"))
    mpi_bcast_cell_structure(CELL_STRUCTURE_NSQUARE);
  else if (ARG1_IS_S("layered")) {
    if (argc > 2) {
      if (!ARG_IS_I(2, n_layers))
	return TCL_ERROR;
      if (n_layers <= 0) {
	Tcl_AppendResult(interp, "layer height should be positive", (char *)NULL);
	return TCL_ERROR;
      }
      determine_n_layers = 0;
    }

    /* check node grid. All we can do is 1x1xn. */
    if (node_grid[0] != 1 || node_grid[1] != 1) {
      node_grid[0] = node_grid[1] = 1;
      node_grid[2] = n_nodes;
      
      err = mpi_bcast_parameter(FIELD_NODEGRID);
    }
    else
      err = 0;

    if (!err)
      mpi_bcast_cell_structure(CELL_STRUCTURE_LAYERED);
  }
  else {
    Tcl_AppendResult(interp, "unkown cell structure type \"", argv[1],"\"", (char *)NULL);
    return TCL_ERROR;
  }
  return gather_runtime_errors(interp, TCL_OK);
}

