/*
  Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project
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
/** \file umbrella_tcl.cpp
 *
 *  Implementation of \ref umbrella_tcl.hpp
 */
#include "umbrella_tcl.hpp"
#include "umbrella.hpp"

#ifdef UMBRELLA

int tclcommand_inter_parse_umbrella(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double k, r;
  int dir;

  if (argc != 4) {
    Tcl_AppendResult(interp, "umbrella needs 3 parameters: "
		     "<k_umbrella> <axis_umbrella> <r0_umbrella>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k)) || (! ARG_IS_I(2, dir)) || (! ARG_IS_D(3, r))) {
    Tcl_AppendResult(interp, "umbrella needs one DOUBLE, one INT, and one DOUBLE parameters: "
		     "<k_umbrella> <axis_umbrella> <r0_umbrella>", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(umbrella_set_params(bond_type, k, dir, r), "bond type must be nonnegative");
}

int tclprint_to_result_umbrellaIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.umbrella.k, buffer);
  Tcl_AppendResult(interp, "UMBRELLA ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.umbrella.dir, buffer);
  Tcl_AppendResult(interp, buffer," ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.umbrella.r, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}

#endif /* UMBREALL */
