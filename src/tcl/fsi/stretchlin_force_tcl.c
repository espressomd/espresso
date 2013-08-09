/*
  Copyright (C) 2010,2011 The ESPResSo project
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

#include "utils.h"
#include "../../tcl/parser.h"
#include "tcl/fsi/stretchlin_force_tcl.h"
#include "fsi/stretchlin_force.h"

/** \file stretchlin_force.h
 *  Routines to calculate the STRETCHLIN_FORCE Energy or/and STRETCHLIN_FORCE force 
 *  for a particle pair. (Dupin2007)
 *  \ref forces.c
*/

/************************************************************/

/// parse parameters for the stretchlin_force potential
int tclcommand_inter_parse_stretchlin_force(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double kslin, r0;

  if (argc != 3) {
    Tcl_AppendResult(interp, "stretchlin_force needs 2 parameters: "
		     "<r0> <kslin>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, r0)) || (! ARG_IS_D(2, kslin)))
    {
      Tcl_AppendResult(interp, "stretchlin_force needs 2 DOUBLE parameters: "
		       "<r0> <kslin>", (char *) NULL);
      return TCL_ERROR;
    }
  
  CHECK_VALUE(stretchlin_force_set_params(bond_type, r0, kslin), "bond type must be nonnegative");
}

int tclprint_to_result_stretchlinforceIA(Tcl_Interp *interp, Bonded_ia_parameters *params){
	char buffer[TCL_DOUBLE_SPACE];
	Tcl_PrintDouble(interp, params->p.stretchlin_force.r0, buffer);
    Tcl_AppendResult(interp, "STRETCHLIN_FORCE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.stretchlin_force.kslin, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
    return (TCL_OK);
}

//#endif

