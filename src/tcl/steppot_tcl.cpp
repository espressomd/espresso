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
/** \file steppot_tcl.c
 *
 *  Implementation of \ref steppot_tcl.hpp
 */
#include "steppot_tcl.hpp"
#include "steppot.hpp"

#ifdef SMOOTH_STEP

int tclprint_to_result_SmStIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_AppendResult(interp, "smooth-step ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_d, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer, "%d", data->SmSt_n);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_eps, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_k0, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  
  return TCL_OK;
}

int tclcommand_inter_parse_SmSt(Tcl_Interp * interp,
				int part_type_a, int part_type_b,
				int argc, char ** argv)
{
  /* parameters needed for LJ */
  double eps, sig, cut, d, k0;
  int n;

  /* get smooth step potential interaction type */
  if (argc < 7) {
    Tcl_AppendResult(interp, "smooth step potential needs 6 parameters: "
		     "<sigma1> <power> <epsilon> <multiplier> <sigma2> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  /* copy smooth step parameters */
  if ((! ARG_IS_D(1, d))     ||
      (! ARG_IS_I(2, n))     ||
      (! ARG_IS_D(3, eps))   ||
      (! ARG_IS_D(4, k0))    ||
      (! ARG_IS_D(5, sig))   ||
      (! ARG_IS_D(6, cut)    )) {
   Tcl_AppendResult(interp, "smooth step potential needs 6 parameters: "
		     "<sigma1> <power> <epsilon> <multiplier> <sigma2> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  if (smooth_step_set_params(part_type_a, part_type_b,
			       d, n, eps, k0, sig,
			       cut) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return 7;
}

#endif
