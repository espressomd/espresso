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
/** \file dihedral_tcl.cpp
 *
 *  Parser for the dihedral potential
 */
#include "dihedral_tcl.hpp"
#include "dihedral.hpp"

/// parse parameters for the dihedral potential
int tclcommand_inter_parse_dihedral(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  int mult;
  double bend, phase;

  if (argc < 4 ) {
    Tcl_AppendResult(interp, "dihedral needs 3 parameters: "
		     "<mult> <bend> <phase>", (char *) NULL);
    return (TCL_ERROR);
  }
  if ( !ARG_IS_I(1, mult) || !ARG_IS_D(2, bend) || !ARG_IS_D(3, phase) ) {
    Tcl_AppendResult(interp, "dihedral needs 3 parameters of types INT DOUBLE DOUBLE: "
		     "<mult> <bend> <phase> ", (char *) NULL);
    return TCL_ERROR;
  }
  
  CHECK_VALUE(dihedral_set_params(bond_type, mult, bend, phase), "bond type must be nonnegative");
}

int tclprint_to_result_dihedralIA(Tcl_Interp *interp,
				  Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  sprintf(buffer, "%d", (int)(params->p.dihedral.mult));
  Tcl_AppendResult(interp, "dihedral ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.dihedral.bend, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.dihedral.phase, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}
