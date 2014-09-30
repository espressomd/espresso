/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "parser.hpp"
#include "interaction_data.hpp"
#include "immersedBoundary/ibm_wall_repulsion.hpp"

/****************
   tclcommand_inter_parse_ibm_wall_repulsion
*****************/

int tclcommand_inter_parse_ibm_wall_repulsion(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  if (!ARG_IS_S_EXACT(1, "cpt") )
  {
    // ****** Setting up a new interaction *********
    
    // Number of arguments
    if ( argc != 2  )
    {
      printf("ibm_wallRep has %d arguments.\n", argc-1);
      Tcl_AppendResult(interp, "Wall repulsion for IBM needs one parameter:\n "
                       "<kappaWall>", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Check correctness of numeric parameters
    double kappaWall;
    if (! ARG_IS_D(1, kappaWall)) {
      Tcl_AppendResult(interp, "kappaWall must be double.", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Setup interaction
    const int status = IBM_WallRepulsion_SetParams(bond_type, kappaWall);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading ibm_wallRep", (char *)NULL);
      return TCL_ERROR;
    }
  }
  else
  {
    // ****** Continue from checkpoint *********
    if ( argc != 2 )
    {
      Tcl_AppendResult(interp, "Wrong number of arguments reading checkpoint for ibm_wallRep ", (char *) NULL);
      return TCL_ERROR;
    }
    // Set interaction
    const int status = IBM_WallRepulsion_ResetParams(bond_type);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading ibm_wallRep checkpoint", (char *)NULL);
      return TCL_ERROR;
    }

  }
}

/***************
   tclprint_to_result_ibm_wall_repulsion
****************/

int tclprint_to_result_ibm_wall_repulsion(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
    // Writes out interaction data to a text file
    // Needed for checkpointing
  
    Tcl_AppendResult(interp, "ibm_wallRep cpt", (char *) NULL);
    return TCL_OK;
}

#endif