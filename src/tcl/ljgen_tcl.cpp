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

/** \file ljgen.c Routines to calculate the generalized lennard jones
 *  energy and/or force for a particle pair. "Generalized" here means
 *  that the LJ energy is of the form
 *
 *  eps * [ b1 * (sigma/(r-r_offset))^a1 - b2 * (sigma/(r-r_offset))^a2 + shift]
 *
 *  \ref forces.c
*/

#include "config.hpp"

#ifdef LENNARD_JONES_GENERIC

#include "ljgen.hpp"
#include "ljgen_tcl.hpp"
#include "communication.hpp"
#include "parser.hpp"

int tclprint_to_result_ljgenIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->LJGEN_eps, buffer);
  Tcl_AppendResult(interp, "lj-gen ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_shift, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  snprintf (buffer, sizeof (buffer), "%d ", data->LJGEN_a1);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  snprintf (buffer, sizeof (buffer), "%d ", data->LJGEN_a2);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  Tcl_PrintDouble(interp, data->LJGEN_b1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_b2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_capradius, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
#ifdef LJGEN_SOFTCORE
  Tcl_PrintDouble(interp, data->LJGEN_lambda, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJGEN_softrad, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
#endif
 
  return TCL_OK;
}


int tclcommand_inter_parse_ljgen(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJGEN */
  double eps, sig, cut, shift, offset, cap_radius, b1, b2;
#ifdef LJGEN_SOFTCORE
  double lambda, softrad;
#endif
  int change, a1, a2;

  /* get lennard-jones interaction type */
  if (argc < 10) {
    Tcl_AppendResult(interp, "lj-gen needs 9 parameters: "
		     "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset> <a1> <a2> <b1> <b2> "
         "[<lj_cap> "
#ifdef LJGEN_SOFTCORE
         "[<lambda> <softrad>]"
#endif
         "]",
		     (char *) NULL);
    return 0;
  }

  /* copy lennard-jones parameters */
  if ((! ARG_IS_D(1, eps))    ||
      (! ARG_IS_D(2, sig))    ||
      (! ARG_IS_D(3, cut))    ||
      (! ARG_IS_D(4, shift))  ||
      (! ARG_IS_D(5, offset)) ||
      (! ARG_IS_I(6, a1))     ||
      (! ARG_IS_I(7, a2))     ||
      (! ARG_IS_D(8, b1))     ||
      (! ARG_IS_D(9, b2))) {
    Tcl_AppendResult(interp, "lj-gen needs 7 DOUBLE and 2 INT parameers: "
		     "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset> <a1> <a2> <b1> <b2> "
         "[<lj_cap> "
#ifdef LJGEN_SOFTCORE
          "[<lambda> <softrad>]"
#endif
          "]",
		     (char *) NULL);
    return ES_ERROR;
  }
  change = 10;
	
  cap_radius = -1.0;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 11 && ARG_IS_D(10, cap_radius))
    change++;
#ifdef LJGEN_SOFTCORE
  lambda = 1.0;
  softrad = 1.0;
  if (argc >= 12 && ARG_IS_D(11, lambda))
    change++;
  else
    Tcl_ResetResult(interp);
  if (argc >= 13 && ARG_IS_D(12, softrad))
    change++;
#endif
  else
    Tcl_ResetResult(interp);
  if (ljgen_set_params(part_type_a, part_type_b,
		       eps, sig, cut, shift, offset, a1, a2, b1, b2,
		       cap_radius
#ifdef LJGEN_SOFTCORE
           , lambda, softrad
#endif
           ) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif
