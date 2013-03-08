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
/** \file subt_lj_tcl.c
 *
 *  Implementation of \ref subt_lj_tcl.h
 */
#include "subt_lj_tcl.h"
#include "subt_lj.h"

#ifdef LENNARD_JONES

int tclcommand_inter_parse_subt_lj(Tcl_Interp *interp, int bond_type,
				   int argc, char **argv)
{
  double k, r;
  if (argc != 3) {
    Tcl_AppendResult(interp, "subt_lj needs 2 dummy parameters: "
		     "<k_subt_lj> <r_subt_lj>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, r))) {
    Tcl_AppendResult(interp, "subt_lj needs 2 dummy DOUBLE parameters: "
		     "<k_subt_lj> <r_subt_lj>", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(subt_lj_set_params(bond_type, k, r), "bond type must be nonnegative");
}

int tclprint_to_result_subt_ljIA(Tcl_Interp *interp,
				 Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.subt_lj.k, buffer);
  Tcl_AppendResult(interp, "SUBT_LJ ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.subt_lj.r, buffer);
  Tcl_AppendResult(interp, buffer,(char *) NULL);
  return TCL_OK;
}

#endif

