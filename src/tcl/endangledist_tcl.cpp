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
/** \file endangledist_tcl.c
 *
 *  Implementation of \ref endangledist_tcl.hpp
 */
#include "utils.hpp"
#include "endangledist.hpp"
#include "endangledist_tcl.hpp"

#ifdef BOND_ENDANGLEDIST

int tclcommand_inter_parse_endangledist(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, phi0, distmin, distmax;

  if (argc != 5) {
    Tcl_AppendResult(interp, "endangledist needs 4 parameters: "
		     "<k> <phi0> <distmin> <distmax>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, bend)) || (! ARG_IS_D(2, phi0)) || (! ARG_IS_D(3, distmin)) || (! ARG_IS_D(4, distmax))) {
    Tcl_AppendResult(interp, "endangledist needs 4 DOUBLE parameters: "
		     "<k> <phi0> <distmin> <distmax> ", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(endangledist_set_params(bond_type, bend, phi0, distmin, distmax), "bond type must be nonnegative");
}

int tclprint_to_result_endangledistIA(Tcl_Interp *interp,
				      Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.endangledist.bend, buffer);
  Tcl_AppendResult(interp, "endangledist ", buffer," ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.endangledist.phi0, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.endangledist.distmin, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.endangledist.distmax, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}

#endif

