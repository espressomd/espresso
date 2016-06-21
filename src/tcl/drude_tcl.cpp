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
/** \file drude_tcl.cpp
 *
 *  Implementation of \ref drude_tcl.hpp
 */

#include "drude_tcl.hpp"
#include "drude.hpp"

#ifdef DRUDE
    
int tclcommand_inter_parse_drude(Tcl_Interp *interp, int bond_type,
				   int argc, char **argv)
{
  double temp_core, gamma_core, temp_drude, gamma_drude, k, mass_drude, r_cut;

  if (argc < 8) {
    Tcl_AppendResult(interp, "drude needs at least 7 parameters: "
		     "<temp_core> <gamma_core> <temp_drude> <gamma_core> <k_drude> <mass_drude> <r_cut> ", (char *) NULL);
    return TCL_ERROR;
  }

  if ( ! ARG_IS_D(1, temp_core) || 
       ! ARG_IS_D(2, gamma_core) ||
       ! ARG_IS_D(3, temp_drude) ||
       ! ARG_IS_D(4, gamma_drude) ||
       ! ARG_IS_D(5, k) ||
       ! ARG_IS_D(6, mass_drude) || 
       ! ARG_IS_D(7, r_cut) ) {
    Tcl_AppendResult(interp, "Parameters should be doubles", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(drude_set_params(bond_type, temp_core, gamma_core, temp_drude, gamma_drude, k, mass_drude, r_cut), "bond type must be nonnegative");
}

int tclprint_to_result_drudeIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.drude.temp_core, buffer);
  Tcl_AppendResult(interp, "DRUDE ", buffer, " ", (char *) NULL);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, params->p.drude.gamma_core, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, params->p.drude.temp_drude, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, params->p.drude.gamma_drude, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, params->p.drude.k, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, params->p.drude.mass_drude, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, params->p.drude.r_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  

  return TCL_OK;
}

#endif

