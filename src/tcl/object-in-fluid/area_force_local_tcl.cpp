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

/** \file area_force_local.h
 *  Routines to calculate the AREA_FORCE_LOCAL energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.c
*/

#include "utils.hpp"
#include "tcl/parser.hpp"
#include "area_force_local_tcl.hpp"
#include "object-in-fluid/area_force_local.hpp"

/************************************************************/


int tclprint_to_result_areaforcelocalIA(Tcl_Interp *interp, Bonded_ia_parameters *params){
	char buffer[TCL_DOUBLE_SPACE];
	Tcl_PrintDouble(interp, params->p.area_force_local.A0_l, buffer);
    Tcl_AppendResult(interp, "AREA_FORCE_LOCAL ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.area_force_local.ka_l, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    return (TCL_OK);
}

/// parse parameters for the area_force_local potential
int tclcommand_inter_parse_area_force_local(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double A0_l, ka_l;

  if (argc != 3) {
    Tcl_AppendResult(interp, "area_force_local needs 2 parameters: "
		     "<A0> <ka>", (char *) NULL);
    return (TCL_ERROR);
  }

 if ((! ARG_IS_D(1, A0_l)) || (! ARG_IS_D(2, ka_l)))
    {
      Tcl_AppendResult(interp, "area_force_local needs 2 DOUBLE parameters: "
		       "<A0_l> <ka_l>", (char *) NULL);
      return TCL_ERROR;
    }

  CHECK_VALUE(area_force_local_set_params(bond_type, A0_l, ka_l), "bond type must be nonnegative");
}

