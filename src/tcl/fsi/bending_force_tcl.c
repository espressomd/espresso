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

/** \file bending_force.h Routines to calculate the bending_force energy or/and
 *  and force for a particle quadruple (two triangles that have 2 particles in common)
*/

#include "../../tcl/parser.h"
#include "fsi/bending_force.h"
#include "tcl/fsi/bending_force_tcl.h"
#include "interaction_data.h"


/// parse parameters for the bending_force potential
int tclcommand_inter_parse_bending_force(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double phi0, kb;

  if (argc != 3 ) {
    Tcl_AppendResult(interp, "bending_force needs 2 parameters: "
		     "<phi0> <kb>", (char *) NULL);
    return (TCL_ERROR);
  }
  if (!ARG_IS_D(1, phi0) || !ARG_IS_D(2, kb) ) {
    Tcl_AppendResult(interp, "bending_force needs 2 parameters of types DOUBLE: "
		     "<phi0> <kb>", (char *) NULL);
    return TCL_ERROR;
  }
  
  CHECK_VALUE(bending_force_set_params(bond_type, phi0, kb), "bond type must be nonnegative");
}

int tclprint_to_result_bendingforceIA(Tcl_Interp *interp, Bonded_ia_parameters *params){
	char buffer[TCL_DOUBLE_SPACE];
	Tcl_PrintDouble(interp, params->p.bending_force.phi0, buffer);
    Tcl_AppendResult(interp, "BENDING_FORCE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.bending_force.kb, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    return (TCL_OK);
}



