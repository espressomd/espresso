/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file quartic_tcl.cpp
 *
 *  Implementation of \ref quartic_tcl.hpp
 */
#include "quartic_tcl.hpp"
#include "quartic.hpp"

int tclcommand_inter_parse_quartic(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double k0, k1, r,r_cut;

  if (argc < 4) {
    Tcl_AppendResult(interp, "quartic needs at least 3 parameters: "
		     "<k0_quartic> <k1_quartic> <r_quartic> [<r_cut>]", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k0)) || (! ARG_IS_D(2, k1)) || (! ARG_IS_D(3, r))) {
    Tcl_AppendResult(interp, "quartic needs at least 3 DOUBLE parameters: "
		     "<k_0c> <k_1> <r_quartic> [<r_cut>]", (char *) NULL);
    return TCL_ERROR;
  }

  if (argc < 5) {
    r_cut = -1.0;
  } else if (! ARG_IS_D(4, r_cut))  {
    Tcl_AppendResult(interp, "<r_cut> should be DOUBLE", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(quartic_set_params(bond_type, k0, k1, r, r_cut), "bond type must be nonnegative");
}

int tclprint_to_result_quarticIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.quartic.k0, buffer);
  Tcl_AppendResult(interp, "QUARTIC ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.quartic.k1, buffer);
  Tcl_AppendResult(interp, buffer," ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.quartic.r, buffer);
  Tcl_AppendResult(interp, buffer," ", (char *) NULL);
  if (params->p.quartic.r_cut > 0.0) {
    Tcl_PrintDouble(interp, params->p.quartic.r_cut, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
  }
  return TCL_OK;
}
