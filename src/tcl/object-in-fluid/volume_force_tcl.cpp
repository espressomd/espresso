/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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

/*  \file volume_force.h
 *  Routines to calculate the VOLUME_FORCE energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.c
*/


#include "utils.hpp"
#include "tcl/parser.hpp"
#include "volume_force_tcl.hpp"
#include "object-in-fluid/volume_force.hpp"

/************************************************************/

int tclprint_to_result_volumeforceIA(Tcl_Interp *interp, Bonded_ia_parameters *params){
	char buffer[TCL_DOUBLE_SPACE];
	Tcl_PrintDouble(interp, params->p.volume_force.V0, buffer);
    Tcl_AppendResult(interp, "VOLUME_FORCE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.volume_force.kv, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    return (TCL_OK);
}

/// parse parameters for the volume_force potential
int tclcommand_inter_parse_volume_force(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double V0, kv;//drvmax;

  if (argc != 3) {
    Tcl_AppendResult(interp, "volume_force needs 2 parameters: "
		     "<V0> <kv>", (char *) NULL);
    return (TCL_ERROR);
  }

 if ((! ARG_IS_D(1, V0)) || (! ARG_IS_D(2, kv)))// || (! ARG_IS_D(3, drvmax)))
    {
      Tcl_AppendResult(interp, "volume needs 2 DOUBLE parameters: "
		       "<V0> <kv>", (char *) NULL);
      return TCL_ERROR;
    }

  CHECK_VALUE(volume_force_set_params(bond_type, V0, kv), "bond type must be nonnegative");
}




