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

/** \file ljcos2.h
 *  Routines to calculate the lennard-jones with cosine tail energy and/or  force 
 *  for a particle pair.  Cosine tail is different from that in ljcos.h
 *  Used for attractive tail/tail interactions in lipid bilayer calculations
 *  \ref forces.c
*/

#include "utils.hpp"

#ifdef LJCOS2
#include <cmath>

#include "ljcos2.hpp"
#include "ljcos2_tcl.hpp"
#include "parser.hpp"
#include "communication.hpp"

int tclprint_to_result_ljcos2IA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->LJCOS2_eps, buffer);
  Tcl_AppendResult(interp, "lj-cos2 ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJCOS2_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJCOS2_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->LJCOS2_w, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  
  return TCL_OK;
}

int tclcommand_inter_parse_ljcos2(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for lj-cos2 */
  double eps, sig, offset, w, cap_radius;
  int change;

  /* get lj-cos2 interaction type */
  if (argc < 5) {
    Tcl_AppendResult(interp, "ljcos2 needs 4 parameters: "
		     "<ljcos2_eps> <ljcos2_sig> <ljcos2_offset> <ljcos2_w>",
		     (char *) NULL);
    return 0;
  }

  /* copy lj-cos2 parameters */
  if ((! ARG_IS_D(1, eps))    ||
      (! ARG_IS_D(2, sig))    ||
      (! ARG_IS_D(3, offset)) ||
      (! ARG_IS_D(4, w))) {
    Tcl_AppendResult(interp, "ljcos2 needs 4 DOUBLE parameters: "
		     "<ljcos2_eps> <ljcos2_sig> <ljcos2_offset> <ljcos2_w>",
		     (char *) NULL);
    return 0;
  }
  change = 5;

  cap_radius = -1;
  /* check wether there is an additional double, cap radius, and parse in */
  if (argc >= 6 && ARG_IS_D(5, cap_radius))
    change++;
  else
    Tcl_ResetResult(interp);

  if (ljcos2_set_params(part_type_a, part_type_b,
			       eps, sig, offset, w
			       ) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif /* ifdef LJCOS2 */
