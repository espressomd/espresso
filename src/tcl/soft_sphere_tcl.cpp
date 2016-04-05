/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file soft_sphere_tcl.cpp
 *
 *  Implementation of \ref soft_sphere_tcl.hpp
 */
#include "soft_sphere_tcl.hpp"
#include "soft_sphere.hpp"

#ifdef SOFT_SPHERE

int tclprint_to_result_softIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->soft_a, buffer);
  Tcl_AppendResult(interp, "soft-sphere ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->soft_n, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->soft_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->soft_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}

int tclcommand_inter_parse_soft(Tcl_Interp * interp,
				int part_type_a, int part_type_b,
				int argc, char ** argv)
{
  /* parameters needed for soft-shere */
  double a, n, cut, offset;
  int change;

  /* get soft-sphere interaction type */
  if (argc < 5) {
    Tcl_AppendResult(interp, "soft-sphere potential needs 4 parameters: "
		     "<soft_a> <soft_n> <soft_cut> <soft_offset>",
		     (char *) NULL);
    return 0;
  }

  /* copy soft-sphere parameters */
  if ((! ARG_IS_D(1, a))     ||
      (! ARG_IS_D(2, n))     ||
      (! ARG_IS_D(3, cut))   ||
      (! ARG_IS_D(4, offset)   )) {
    Tcl_AppendResult(interp, "soft-sphere potential needs 4 parameters: "
		     "<soft_a> <soft_n> <soft_cut> <soft_offset>",
		     (char *) NULL);
    return 0;
  }
  change = 5;
	
  
  Tcl_ResetResult(interp);
  if (soft_sphere_set_params(part_type_a, part_type_b,
                             a, n, cut, offset) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

#endif

