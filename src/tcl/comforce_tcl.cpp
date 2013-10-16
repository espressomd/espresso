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
/** \file comforce_tcl.cpp
 *
 *  Implementation of \ref comforce_tcl.hpp.
 */
#include "utils.hpp"
#include "parser.hpp"
#include "comforce_tcl.hpp"
#include "comforce.hpp"

#ifdef COMFORCE

int tclprint_to_result_comforceIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%d",data->COMFORCE_flag);
  Tcl_AppendResult(interp, "comforce ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",data->COMFORCE_dir);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->COMFORCE_force, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->COMFORCE_fratio, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    
  return TCL_OK;
}

int tclcommand_inter_parse_comforce(Tcl_Interp * interp,
				    int part_type_a, int part_type_b,
				    int argc, char ** argv)
{
  int flag, dir, change; 
  double force, fratio;

  if (argc != 5) {
    Tcl_AppendResult(interp, "comforce needs 4 parameters: "
		     "<comforce_flag> <comforce_dir> <comforce_force> <comforce_fratio>",
		     (char *) NULL);
    return 0;
  }
  
  if (part_type_a == part_type_b) {
    Tcl_AppendResult(interp, "comforce needs 2 different types ", (char *) NULL);
    return 0;
  }

  /* copy comforce parameters */
  if ((! ARG_IS_I(1, flag)) || (! ARG_IS_I(2, dir)) || (! ARG_IS_D(3, force)) || (! ARG_IS_D(4, fratio)) ) {
    Tcl_AppendResult(interp, "comforce needs 2 INTEGER 1 DOUBLE parameter: "
		     "<comforce_flag> <comforce_dir> <comforce_force> <comforce_fratio>", (char *) NULL);
    return 0;
  }
    
  change = 5;
    
  switch (comforce_set_params(part_type_a, part_type_b, flag, dir, force, fratio)) {
  case 1:
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  case 2:
    Tcl_AppendResult(interp, "works only with a single CPU", (char *) NULL);
    return 0;
  }
  return change;
}

#endif
