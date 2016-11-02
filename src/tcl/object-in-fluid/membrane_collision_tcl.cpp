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
/** \file membrane_collision_tcl.cpp
 *
 *  Implementation of \ref membrane_collision_tcl.hpp
 */
#include "membrane_collision_tcl.hpp"
#include "../../core/object-in-fluid/membrane_collision.hpp"

#ifdef MEMBRANE_COLLISION

int tclprint_to_result_membraneIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_AppendResult(interp, "membrane ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->membrane_a, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->membrane_n, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->membrane_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->membrane_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}

int tclcommand_inter_parse_membrane(Tcl_Interp * interp,
				int part_type_a, int part_type_b,
				int argc, char ** argv)
{
  /* parameters needed for membrane collision */
  double a, n, offset, cut;
  int change;

  if (argc < 5) {
    Tcl_AppendResult(interp, "membrane collision potential needs 4 parameters: "
		     "<membrane_a> <membrane_n> <membrane_cut> <membrane_offset>",
		     (char *) NULL);
    return 0;
  }

  /* copy membrane collision parameters */
    if ((! ARG_IS_D(1, a))     ||
        (! ARG_IS_D(2, n))     ||
        (! ARG_IS_D(3, cut))   ||
        (! ARG_IS_D(4, offset)   )) {
        Tcl_AppendResult(interp, "membrane collision potential needs 4 parameters: "
                         "<soft_a> <soft_n> <soft_cut> <soft_offset>",
		     (char *) NULL);
    return 0;
  }
  change = 5;

  if (membrane_collision_set_params(part_type_a, part_type_b,
                             a, n, cut, offset) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif

