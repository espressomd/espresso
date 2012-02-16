/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file mmm1d.c  MMM1D algorithm for long range coulomb interaction.
 *
 *  For more information about MMM1D, see \ref mmm1d.h "mmm1d.h".
 */

#include <mpi.h>
#include <tcl.h>
#include "utils.h"
#include "mmm1d.h"
#include "polynom.h"
#include "specfunc.h"
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "tuning.h"
#include "interaction_data.h"
#include "mmm-common.h"
#include "parser.h"
#include "errorhandling.h"

#ifdef ELECTROSTATICS

int tclprint_to_result_MMM1D(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, sqrt(mmm1d_params.far_switch_radius_2), buffer);
  Tcl_AppendResult(interp, "mmm1d ", buffer, " ",(char *) NULL);
  sprintf(buffer, "%d", mmm1d_params.bessel_cutoff);
  Tcl_AppendResult(interp, buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, mmm1d_params.maxPWerror, buffer);
  Tcl_AppendResult(interp, buffer,(char *) NULL);

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_mmm1d(Tcl_Interp *interp, int argc, char **argv)
{
  double switch_rad, maxPWerror;
  int bessel_cutoff;

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb mmm1d <switch radius> "
		     "{<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
    return TCL_ERROR;
  }

  if (ARG0_IS_S("tune")) {
    /* autodetermine bessel cutoff AND switching radius */
    if (! ARG_IS_D(1, maxPWerror))
      return TCL_ERROR;
    bessel_cutoff = -1;
    switch_rad = -1;
  }
  else {
    if (argc == 2) {
      /* autodetermine bessel cutoff */
      if ((! ARG_IS_D(0, switch_rad)) ||
	  (! ARG_IS_D(1, maxPWerror))) 
	return TCL_ERROR;
      bessel_cutoff = -1;
    }
    else if (argc == 3) {
      /* fully manual */
      if((! ARG_IS_D(0, switch_rad)) ||
	 (! ARG_IS_I(1, bessel_cutoff)) ||
	 (! ARG_IS_D(2, maxPWerror))) 
	return TCL_ERROR;

      if (bessel_cutoff <=0) {
	Tcl_AppendResult(interp, "bessel cutoff too small", (char *)NULL);
	return TCL_ERROR;
      }
    }
    else {
      Tcl_AppendResult(interp, "wrong # arguments: inter coulomb mmm1d <switch radius> "
		       "{<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
      return TCL_ERROR;
    }
    
    if (switch_rad <= 0 || switch_rad > box_l[2]) {
      Tcl_AppendResult(interp, "switching radius is not between 0 and box_l[2]", (char *)NULL);
      return TCL_ERROR;
    }
  }

  MMM1D_set_params(switch_rad, bessel_cutoff, maxPWerror);

  char *log = NULL;
  int result = mmm1d_tune(&log) == ES_OK ? TCL_OK : TCL_ERROR;

  Tcl_AppendResult(interp, log, NULL);
  if (log)
    free(log);

  return gather_runtime_errors(interp, result);
}

#endif
