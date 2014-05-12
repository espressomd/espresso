/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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

/** \file morse_tcl.cpp
 *  TCL interface for the Morse potential
*/

#include "morse_tcl.hpp"

#ifdef MORSE

#include "morse.hpp"
#include "interaction_data.hpp"
#include "parser.hpp"
#include "forcecap_tcl.hpp"

int tclprint_to_result_morseIA(Tcl_Interp *interp, int i, int j)
{
    char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->MORSE_eps, buffer);
  Tcl_AppendResult(interp, "morse ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->MORSE_alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->MORSE_rmin, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->MORSE_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->MORSE_capradius, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

/// parser for the forcecap
int tclcommand_inter_parse_morseforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  if (argc==1) {
    fprintf(stderr, "WARNING: \"inter morseforcecap\" is deprecated "
	    "and will be removed in some further version. "
	    "Use \"inter forcecap\" instead.\n");
  }
  return tclcommand_inter_parse_forcecap(interp, argc, argv);
}


int tclcommand_inter_parse_morse(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for MORSE */
  double eps, alpha, rmin, cut, cap_radius;
  int change;

  /* get morse interaction type */
  if (argc < 5) {
    Tcl_AppendResult(interp, "morse needs 4 parameters: "
		     "<morse_eps> <morse_alpha> <morse_rmin> <morse_cut>",
		     (char *) NULL);
    return 0;
  }

  /* copy morse parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, alpha))   ||
      (! ARG_IS_D(3, rmin))   ||
      (! ARG_IS_D(4, cut)   )) {
    Tcl_AppendResult(interp, "morse needs 4 DOUBLE parameters: "
		     "<morse_eps> <morse_alpha> <morse_rmin> <morse_cut>",
		     (char *) NULL);
    return 0;
  }
  change = 5;
	
  cap_radius = -1.0;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 6 && ARG_IS_D(5, cap_radius))
    change++;
  else
    Tcl_ResetResult(interp);
  if (morse_set_params(part_type_a, part_type_b,
			       eps, alpha, rmin, cut, cap_radius) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif
