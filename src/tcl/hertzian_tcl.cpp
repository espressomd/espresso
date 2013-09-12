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
/** \file hertzian_tcl.c
 *
 *  Implementation of \ref hertzian_tcl.hpp
 */
#include "hertzian_tcl.hpp"
#include "hertzian.hpp"

#ifdef HERTZIAN

int tclprint_to_result_HertzianIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->Hertzian_eps, buffer);
  Tcl_AppendResult(interp, "hertzian ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->Hertzian_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_parse_hertzian(Tcl_Interp * interp,
				    int part_type_a, int part_type_b,
				    int argc, char ** argv)
{
  /* parameters needed for Hertzian */
  double eps, sig;

  /* get interaction type */
  if (argc < 3) {
    Tcl_AppendResult(interp, "Hertzian potential needs 2 parameters: "
		     "<epsilon> <sigma>", (char *) NULL);
    return 0;
  }

  /* copy parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, sig))) {
    Tcl_AppendResult(interp, "Hertzian potential needs 2 parameters: "
		     "<epsilon> <sigma>", (char *) NULL);
    return 0;
  }
  
  Tcl_ResetResult(interp);
  if (hertzian_set_params(part_type_a, part_type_b,
			  eps, sig) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  /* return number of used parameters */
  return 3;
}

#endif
