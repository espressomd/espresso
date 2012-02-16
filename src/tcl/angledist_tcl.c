/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file angledist_tcl.c
 *
 *  Implementation of \ref angledist_tcl.h
 */
#include "angledist_tcl.h"

#ifdef BOND_ANGLEDIST
#include "angledist.h"

int tclcommand_inter_parse_angledist(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, phimin, distmin, phimax, distmax;

  if (argc != 6) {
    Tcl_AppendResult(interp, "angledist needs 5 parameters: "
		     "<bend> <phimin> <distmin> <phimax> <distmax>", (char *) NULL);
    printf ("argc=%d\n",argc);
    return (TCL_ERROR);
  }

  if (! ARG_IS_D(1, bend)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<bend> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(2, phimin)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
                     "<phimin> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(3, distmin)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<distmin> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(4, phimax)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
                     "<phimax> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(5, distmax)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<distmax> ", (char *) NULL);
    return TCL_ERROR;
  }


  CHECK_VALUE(angledist_set_params(bond_type, bend, phimin, distmin, phimax, distmax), "bond type must be nonnegative");
}

int tclprint_to_result_angledistIA(Tcl_Interp *interp,
				   Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.angledist.bend, buffer);
  Tcl_AppendResult(interp, "angledist ", buffer," ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.angledist.phimin, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.angledist.distmin, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.angledist.phimax, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.angledist.distmax, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}

#endif
