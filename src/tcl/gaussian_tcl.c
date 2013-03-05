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
/** \file gaussian_tcl.c
 *
 *  Implementation of \ref gaussian_tcl.h
 */
#include "gaussian_tcl.h"
#include "gaussian.h"

#ifdef GAUSSIAN

int tclprint_to_result_GaussianIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->Gaussian_eps, buffer);
  Tcl_AppendResult(interp, "gaussian ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->Gaussian_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->Gaussian_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_parse_gaussian(Tcl_Interp * interp,
				    int part_type_a, int part_type_b,
				    int argc, char ** argv)
{
  /* parameters needed for Gaussian */
  double eps, sig, cut;

  /* copy parameters */
  if ((argc < 4) ||
      (! ARG_IS_D(1, eps)) ||
      (! ARG_IS_D(2, sig)) || 
      (! ARG_IS_D(3, cut))) {
    Tcl_AppendResult(interp, "Gaussian potential needs 3 parameters: "
		     "<epsilon> <sigma> <cut>", (char *) NULL);
    return 0;
  }
  
  Tcl_ResetResult(interp);
  if (gaussian_set_params(part_type_a, part_type_b,
			  eps, sig, cut) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  /* return number of used parameters */
  return 4;
}

#endif
