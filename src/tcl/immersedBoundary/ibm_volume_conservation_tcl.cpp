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
#include "immersedBoundary/ibm_volume_conservation_tcl.hpp"
#include "immersedBoundary/ibm_volume_conservation.hpp"

/****************
   tclcommand_inter_parse_ibm_volcons
*****************/

int tclcommand_inter_parse_ibm_volume_conservation(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  if (!ARG_IS_S_EXACT(1, "cpt") )
  {
    // ****** Setting up a new interaction *********
    
    // Number of arguments
    if ( argc != 3  )
    {
      printf("ibm_volcons has %d arguments.\n", argc-1);
      Tcl_AppendResult(interp, "VolCons for IBM needs two parameters:\n "
                       "<softID> <kappaV>", (char *) NULL);
      return TCL_ERROR;
    }
    
    int softID;
    double kappaV;
    
    // Check correctness of numeric parameters
    if (! ARG_IS_I(1, softID)) {
      Tcl_AppendResult(interp, "softID must be int\n", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (! ARG_IS_D(2, kappaV)) {
      Tcl_AppendResult(interp, "kappaV must be double\n", (char *) NULL);
      return TCL_ERROR;
    }
    
/*    bool writeCom = false;
    if ( argc == 4 )
    {
      if ( ARG_IS_S_EXACT(3, "writeCOM") ) writeCom = true;
      else
      {
        Tcl_AppendResult(interp, "illegal optional argument: must be writeCOM or nothing", (char *) NULL);
        return TCL_ERROR;
      }
    }*/
    
    // Setup interaction
    const int status = IBM_VolumeConservation_SetParams(bond_type, softID, kappaV);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading volcons", (char *)NULL);
      return TCL_ERROR;
    }
  }
  else
  {
    // ****** Continue from checkpoint *********
    if ( argc != 2 )
    {
      Tcl_AppendResult(interp, "Wrong number of arguments reading checkpoint for ibm_volcons ", (char *) NULL);
      return TCL_ERROR;
    }
    
    // ****** Continue from checkpoint *********
    double volRef;
    if (! ARG_IS_D(2, volRef))
    {
      Tcl_AppendResult(interp, "volRef must be double ", (char *) NULL);
      return TCL_ERROR;
    }
    
    // Reset interaction
    const int status = IBM_VolumeConservation_ResetParams(bond_type, volRef);
    if (status == ES_ERROR)
    {
      Tcl_AppendResult(interp, "Unspecified error reading volcons checkpoint", (char *)NULL);
      return TCL_ERROR;
    }
  }
}

/***************
   tclprint_to_result_ibm_volcons
****************/
    
int tclprint_to_result_ibm_volume_conservation(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  // Writes out interaction data to a text file
  // Needed for checkpointing
  
  char buffer[TCL_DOUBLE_SPACE];
  
  Tcl_PrintDouble(interp, params->p.ibmVolConsParameters.volRef, buffer);
  Tcl_AppendResult(interp, "ibm_volcons checkpoint ", buffer," ", (char *) NULL);
  
  return TCL_OK;
}

#endif