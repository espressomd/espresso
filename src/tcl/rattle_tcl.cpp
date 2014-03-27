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
/** \file rattle_tcl.cpp
 *
 *  Implementation of \ref rattle_tcl.hpp
 */
#include "parser.hpp"
#include "integrate.hpp"
#include "rattle.hpp"
#include "interaction_data.hpp"

#ifdef BOND_CONSTRAINT

/// parse parameters for the rigid bonds
int tclcommand_inter_parse_rigid_bond(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double d, p_tol, v_tol;
    
  if (argc != 4) {
    Tcl_AppendResult(interp, "rigid bond needs 3 parameters: "
		     "<constrained_bond_distance> <Positional_tolerance> <Velocity_tolerance>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, d)) || (! ARG_IS_D(2, p_tol)) || (! ARG_IS_D(3, v_tol)) ) {
    Tcl_AppendResult(interp, "rigid bond needs 3 DOUBLE parameters: "
		     "<constrained_bond_distance> <Positional_tolerance> <Velocity_tolerance>", (char *) NULL);
    return TCL_ERROR;
  }

  if (time_step < 0. ) {
    Tcl_AppendResult(interp, "set time_step before declaring rigid_bond" , (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(rigid_bond_set_params(bond_type, d, p_tol, v_tol), "bond type must be nonnegative");
}

int tclprint_to_result_rigid_bond(Tcl_Interp *interp,
				  Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, sqrt(params->p.rigid_bond.d2), buffer);
  Tcl_AppendResult(interp, "RIGID_BOND ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.rigid_bond.p_tol/2.0, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.rigid_bond.v_tol/time_step, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}

#endif
