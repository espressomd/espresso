/*
  Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project
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

/** \file cos2.hpp
 *  Routines to calculate a flat potential with cosine tail energy and/or  force 
 *  for a particle pair.  Cosine tail is different from that in ljcos.hpp
 *  Used for attractive tail/tail interactions in lipid bilayer calculations.
 *  Same potential as ljcos2 without Lennard-Jones part.
 *  \ref forces.cpp
*/

#include "utils.hpp"

#ifdef COS2
#include <cmath>

#include "cos2.hpp"
#include "cos2_tcl.hpp"
#include "parser.hpp"
#include "communication.hpp"

int tclprint_to_result_cos2IA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->COS2_eps, buffer);
  Tcl_AppendResult(interp, "cos2 ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->COS2_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->COS2_w, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  
  return TCL_OK;
}

int tclcommand_inter_parse_cos2(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for cos2 */
  double eps, offset, w;
  int change;

  /* get cos2 interaction type */
  if (argc < 4) {
    Tcl_AppendResult(interp, "cos2 needs 3 parameters: "
		     "<cos2_eps> <cos2_offset> <cos2_w>",
		     (char *) NULL);
    return 0;
  }

  /* copy cos2 parameters */
  if ((! ARG_IS_D(1, eps))    ||
      (! ARG_IS_D(2, offset)) ||
      (! ARG_IS_D(3, w))) {
    Tcl_AppendResult(interp, "cos2 needs 3 DOUBLE parameters: "
		     "<cos2_eps> <cos2_offset> <cos2_w>",
		     (char *) NULL);
    return 0;
  }
  change = 4;

  if (cos2_set_params(part_type_a, part_type_b,
			       eps, offset, w
			       ) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif /* ifdef COS2 */
