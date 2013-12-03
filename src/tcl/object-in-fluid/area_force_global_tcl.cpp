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

/** \file area_force_global.hpp
 *  Routines to calculate the AREA_FORCE_GLOBAL energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "tcl/parser.hpp"
#include "area_force_global_tcl.hpp"
#include "object-in-fluid/area_force_global.hpp"

int tclprint_to_result_areaforceglobalIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, params->p.area_force_global.A0_g, buffer);
    Tcl_AppendResult(interp, "AREA_FORCE_GLOBAL ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.area_force_global.ka_g, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  return (TCL_OK);
}

/// parse parameters for the area_force_global potential
int tclcommand_inter_parse_area_force_global(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double A0_g, ka_g;

  if (argc != 3) {
    Tcl_AppendResult(interp, "area_force_global needs 2 parameters: "
		     "<A0_g> <ka_g>", (char *) NULL);
    return (TCL_ERROR);
  }

 if ((! ARG_IS_D(1, A0_g)) || (! ARG_IS_D(2, ka_g)))
    {
      Tcl_AppendResult(interp, "area_force_global needs 2 DOUBLE parameters: "
		       "<A0_g> <ka_g>", (char *) NULL);
      return TCL_ERROR;
    }

  CHECK_VALUE(area_force_global_set_params(bond_type, A0_g, ka_g), "bond type must be nonnegative");
}

