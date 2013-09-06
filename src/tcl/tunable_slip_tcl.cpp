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
/** \file tunable_slip_tcl.c
 *
 *  Implementation of \ref tunable_slip_tcl.hpp
 */
#include "parser.hpp"

#ifdef TUNABLE_SLIP
#include "tunable_slip.hpp"
#include "interaction_data.hpp"

int tclprint_to_result_tunable_slipIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_temp, buffer);
  Tcl_AppendResult(interp, "tunable_slip ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_gamma, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_r_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_time, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_vx, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_vy, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->TUNABLE_SLIP_vz, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}

int tclcommand_inter_parse_tunable_slip(Tcl_Interp * interp,
					int part_type_a, int part_type_b,
					int argc, char ** argv)
{
  double temp, gamma, r_cut, time, vx, vy, vz;
  int change;

  if (argc < 8) {
    Tcl_AppendResult(interp, "tunable_slip needs 7 parameters: "
		     "<tunable_slip_temp> <tunable_slip_gamma> <tunable_slip_r_cut> <tunable_slip_time> <vx> <vy> <vz>",
		     (char *) NULL);
    return 0;
  }

  /* copy lj-cos parameters */
  if ((! ARG_IS_D(1, temp))    ||
      (! ARG_IS_D(2, gamma))   ||
      (! ARG_IS_D(3, r_cut))   ||
      (! ARG_IS_D(4, time))   ||
      (! ARG_IS_D(5, vx))   ||
      (! ARG_IS_D(6, vy))   ||
      (! ARG_IS_D(7, vz)    )) {
    Tcl_AppendResult(interp, "tunable_slip needs 7 DOUBLE parameters: "
		     "<tunable_slip_temp> <tunable_slip_gamma> <tunable_slip_r_cut> <tunable_slip_time> <vx> <vy> <vz>",
		     (char *) NULL);
    return 0;
  }
  change = 8;

  if (tunable_slip_set_params(part_type_a, part_type_b, temp, gamma, r_cut, time, vx, vy, vz) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }

  return change;
}

#endif
