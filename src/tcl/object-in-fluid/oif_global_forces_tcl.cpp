/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
  
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

/** \file oif_global_forces.hpp
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "tcl/parser.hpp"
#include "oif_global_forces_tcl.hpp"
#include "object-in-fluid/oif_global_forces.hpp"

int tclprint_to_result_oifglobalforcesIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, params->p.oif_global_forces.A0_g, buffer);
    Tcl_AppendResult(interp, "OIF_GLOBAL_FORCES ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_global_forces.ka_g, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_global_forces.V0, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, params->p.oif_global_forces.kv, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  return (TCL_OK);
}

/// parse parameters for the oif_global_forces potential
int tclcommand_inter_parse_oif_global_forces(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double A0_g, ka_g, V0, kv;

  if (argc != 5) {
    Tcl_AppendResult(interp, "oif_global_forces needs 4 parameters: "
		     "<A0_g> <ka_g> <V0> <kv>", (char *) NULL);
    return (TCL_ERROR);
  }

 if ((! ARG_IS_D(1, A0_g)) || (! ARG_IS_D(2, ka_g)) || (! ARG_IS_D(3, V0)) || (! ARG_IS_D(4, kv)))
    {
      Tcl_AppendResult(interp, "oif_global_forces needs 4 DOUBLE parameters: "
		       "<A0_g> <ka_g>  <V0> <kv>", (char *) NULL);
      return TCL_ERROR;
    }

  CHECK_VALUE(oif_global_forces_set_params(bond_type, A0_g, ka_g, V0, kv), "bond type must be nonnegative");
}

