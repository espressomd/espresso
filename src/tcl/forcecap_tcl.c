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
#include "utils.h"
#include "lj.h"
#include "parser.h"
#include "mol_cut.h"
#include "communication.h"
#include "forcecap_tcl.h"

/// parser for the forcecap
int tclcommand_inter_parse_forcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE];

  double forcecap;
  
  if (argc == 0) {
    if (force_cap == -1.0)
      Tcl_AppendResult(interp, "forcecap individual", (char *) NULL);
    else {
      Tcl_PrintDouble(interp, force_cap, buffer);
      Tcl_AppendResult(interp, "forcecap ", buffer, (char *) NULL);
    }
    return TCL_OK;
  }

  if (argc > 1) {
    Tcl_AppendResult(interp, "inter forcecap takes at most 1 parameter",
		     (char *) NULL);      
    return TCL_ERROR;
  }
  
  if (ARG0_IS_S("individual")){
    forcecap = -1.0;
    CHECK_VALUE(forcecap_set_params(forcecap),
	      "If you can read this, you should change it. (Use the source Luke!)");
   }
  else if (! ARG0_IS_D(forcecap) || forcecap < 0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(forcecap_set_params(forcecap),
	      "If you can read this, you should change it. (Use the source Luke!)");
}
