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
/** \file gb_tcl.c
 *
 *  Implementation of \ref gb_tcl.h
 */
#include "gb_tcl.h"

#ifdef GAY_BERNE
#include "gb.h"
#include "communication.h"

int tclprint_to_result_gbIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->GB_eps, buffer);
  Tcl_AppendResult(interp, "gay-berne ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_k1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_k2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_mu, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_nu, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_chi1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->GB_chi2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_parse_gb(Tcl_Interp * interp,
			      int part_type_a, int part_type_b,
			      int argc, char ** argv)
{
  double tmp;
  double eps, sig, cut;
  double k1, k2, mu, nu;
  int change;

  /* there are 9 parameters for gay-berne, but you read in only 7 of them.
     The rest is calculated in gay_berne_set_params.
  */

  if (argc < 8) {
    Tcl_AppendResult(interp, "gay-berne needs 7 parameters: "
		     "<gb_eps> <gb_sig> <gb_cut> <gb_k1> <gb_k2> <gb_mu> <gb_nu>",
		     (char *) NULL);
    return 0;
  }

  /* copy gay-berne parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, sig))   ||
      (! ARG_IS_D(3, cut))   ||
      (! ARG_IS_D(4, k1 ))   ||
      (! ARG_IS_D(5, k2 ))   ||
      (! ARG_IS_D(6, mu ))   ||	
      (! ARG_IS_D(7, nu )    )) {
    Tcl_AppendResult(interp, "gay-berne needs 7 DOUBLE parameters: "
		     "<gb_eps> <gb_sig> <gb_cut> <gb_k1> <gb_k2> <gb_mu> <gb_nu>",
		     (char *) NULL);
    return 0;
  }
  change = 8;

  if (argc >= 10 && ARG_IS_D(8, tmp) && ARG_IS_D(9, tmp))
    change += 2;
  else
    Tcl_ResetResult(interp);

  if (gay_berne_set_params(part_type_a, part_type_b, eps, sig, cut, k1, k2, mu, nu) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif
