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
/** \file angle_cossquare_tcl.c
 *
 *  Implementation of \ref angle_cossquare_tcl.hpp
 */
#include "angle_cossquare_tcl.hpp"

#ifdef BOND_ANGLE
#include "angle_cossquare.hpp"

#include "communication.hpp"

/// parse parameters for the angle_cossquare potential
int tclcommand_inter_parse_angle_cossquare(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, phi0;

  /* the optional parameter phi0 is due to backwards compatibility and is set to PI if not given */
  if (argc != 2 && argc != 3) {
    Tcl_AppendResult(interp, "angle_cossquare needs 1 or 2 parameters: "
		     "<bend> [<phi0>]", (char *) NULL);
    return (TCL_ERROR);
  }

  if (! ARG_IS_D(1, bend)) {
    Tcl_AppendResult(interp, "angle_cossquare needs a DOUBLE parameter: "
		     "<bend> ", (char *) NULL);
    return TCL_ERROR;
  }

  /* special treatment of the optional parameter phi0 */
  if (argc == 3) {
    if (! ARG_IS_D(2, phi0)) {
      Tcl_AppendResult(interp, "angle_cossquare needs a DOUBLE parameter: "
		       "<phi0> ", (char *) NULL);
      return TCL_ERROR;
    }
  } else {
    phi0 = PI;
  }
  CHECK_VALUE(angle_cossquare_set_params(bond_type, bend, phi0), "bond type must be nonnegative");
}

int tclprint_to_result_angle_cossquareIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, params->p.angle_cossquare.bend, buffer);
  Tcl_AppendResult(interp, "angle_cossquare ", buffer," ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.angle_cossquare.phi0, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}

#endif

