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
/** \file mmm1d.cpp  MMM1D algorithm for long range coulomb interaction.
 *
 *  For more information about MMM1D, see \ref mmm1d.hpp "mmm1d.h".
 */

#include <mpi.h>
#include <tcl.h>
#include "utils.hpp"
#include "mmm1d.hpp"
#include "polynom.hpp"
#include "specfunc.hpp"
#include "communication.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "tuning.hpp"
#include "interaction_data.hpp"
#include "mmm-common.hpp"
#include "parser.hpp"
#include "errorhandling.hpp"

#ifdef ELECTROSTATICS

int tclprint_to_result_MMM1D(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, sqrt(mmm1d_params.far_switch_radius_2), buffer);
  Tcl_AppendResult(interp, "mmm1d ", buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, mmm1d_params.maxPWerror, buffer);
  Tcl_AppendResult(interp, buffer,(char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_mmm1d(Tcl_Interp *interp, int argc, char **argv)
{
  double switch_rad, maxPWerror;

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb mmm1d <switch radius> "
		     "{<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
    return TCL_ERROR;
  }

  if (ARG0_IS_S("tune")) {
    /* autodetermine bessel cutoff AND switching radius */
    if (! ARG_IS_D(1, maxPWerror))
      return TCL_ERROR;
    switch_rad = -1;
  }
  else {
    if (argc == 2) {
      /* autodetermine bessel cutoff */
      if ((! ARG_IS_D(0, switch_rad)) ||
	  (! ARG_IS_D(1, maxPWerror))) 
	return TCL_ERROR;
    }
    else {
      Tcl_AppendResult(interp, "wrong # arguments: inter coulomb mmm1d <switch radius> "
		       "<maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (switch_rad <= 0 || switch_rad > box_l[2]) {
      Tcl_AppendResult(interp, "switching radius is not between 0 and box_l[2]", (char *)NULL);
      return TCL_ERROR;
    }
  }

  MMM1D_set_params(switch_rad, maxPWerror);

  char *log = NULL;
  int result = mmm1d_tune(&log) == ES_OK ? TCL_OK : TCL_ERROR;

  Tcl_AppendResult(interp, log, NULL);
  if (log)
    free(log);

  return gather_runtime_errors(interp, result);
}

#endif
