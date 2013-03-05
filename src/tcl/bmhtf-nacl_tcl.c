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
/** \file bmhtf-nacl_tcl.c
 *
 *  Implementation of \ref bmhtf-nacl_tcl.h
 */
#include "bmhtf-nacl_tcl.h"

#ifdef BMHTF_NACL
#include "bmhtf-nacl.h"
#include "communication.h"

int tclprint_to_result_BMHTFIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_AppendResult(interp, "bmhtf-nacl ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_A, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_B, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_C, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_D, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  
  return TCL_OK;
}

int tclcommand_inter_parse_BMHTF(Tcl_Interp * interp,
				 int part_type_a, int part_type_b,
				 int argc, char ** argv)
{
  double A, B, C, D, sig, cut;

  if (argc < 7) {
    Tcl_AppendResult(interp, "BMHTF NaCl potential needs 6 parameters: "
		     "<A> <B> <C> <D> <sigma> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  /* copy smooth step parameters */
  if ((! ARG_IS_D(1, A))    ||
      (! ARG_IS_D(2, B))    ||
      (! ARG_IS_D(3, C))    ||
      (! ARG_IS_D(4, D))    ||
      (! ARG_IS_D(5, sig))  ||
      (! ARG_IS_D(6, cut)   )) {
    Tcl_AppendResult(interp, "BMHTF NaCl potential needs 6 parameters: "
		     "<A> <B> <C> <D> <sigma> <cutoff>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (BMHTF_set_params(part_type_a, part_type_b,
		       A, B, C, D, sig, cut) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return 7;
}

#endif
