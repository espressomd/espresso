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
/** \file bonded_coulomb_tcl.cpp
 *
 *  Implementation of \ref bonded_coulomb_tcl.hpp
 */
#include "bonded_coulomb_tcl.hpp"
#include "bonded_coulomb.hpp"

#ifdef ELECTROSTATICS

int tclcommand_inter_parse_bonded_coulomb(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double prefactor;

  if (argc != 2) {
    Tcl_AppendResult(interp, "bonded_coulomb exactly one parameter: "
		     "<prefactor>", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(1, prefactor)) {
    Tcl_AppendResult(interp, "bonded_coulomb needs 1 DOUBLE parameters: "
		     "<prefactor>", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(bonded_coulomb_set_params(bond_type, prefactor), "bond type must be nonnegative");
}

int tclprint_to_result_bonded_coulombIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, params->p.bonded_coulomb.prefactor, buffer);
  Tcl_AppendResult(interp, "BONDED_COULOMB ", buffer, " ", (char *) NULL);

  return TCL_OK;
}

#endif
