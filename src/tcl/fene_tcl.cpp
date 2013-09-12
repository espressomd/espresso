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
/** \file fene_tcl.c
 *
 *  Implementation of \ref fene_tcl.hpp
 */
#include "utils.hpp"
#include "parser.hpp"
#include "fene_tcl.hpp"
#include "fene.hpp"

int tclprint_to_result_feneIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, params->p.fene.k, buffer);
  Tcl_AppendResult(interp, "FENE ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.fene.drmax, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  Tcl_PrintDouble(interp, params->p.fene.r0, buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  return (TCL_OK);
}

int tclcommand_inter_parse_fene(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double k, drmax, r0;

  if (argc != 3 && argc != 4) {
    Tcl_AppendResult(interp, "fene needs 2 or 3 parameters: "
		     "<k> <drmax> [<r0>]", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, drmax)))
    {
      Tcl_AppendResult(interp, "fene needs 2 or 3 DOUBLE parameters: "
		       "<k> <drmax> [<r0>]", (char *) NULL);
      return TCL_ERROR;
    }

  if (argc == 4) {
    if (! ARG_IS_D(3, r0))
      {
	Tcl_AppendResult(interp, "fene needs 2 or 3 DOUBLE parameters: "
			 "<k> <drmax> [<r0>]", (char *) NULL);
	return TCL_ERROR;
      }
  } else {
    /* default value for r0 is 0.0. */
    r0 = 0.0;
  }
  
  CHECK_VALUE(fene_set_params(bond_type, k, drmax, r0), "bond type must be nonnegative");
}

