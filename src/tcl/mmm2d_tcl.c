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
/** \file mmm2d.c  MMM2D algorithm for long range coulomb interaction.
 *
 *  For more information about MMM2D, see \ref mmm2d.h "mmm2d.h".
 */

//#include <math.h>
//#include <mpi.h>
#include "utils.h"
#include "communication.h"
//#include "particle_data.h"
//#include "interaction_data.h"
#include "cells.h"
#include "mmm2d.h"
//#include "mmm-common.h"
//#include "specfunc.h"
//#include "integrate.h"
//#include "layered.h"
#include "parser.h"

#ifdef ELECTROSTATICS

int tclprint_to_result_MMM2D(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, mmm2d_params.maxPWerror, buffer);
  Tcl_AppendResult(interp, "mmm2d ", buffer,(char *) NULL);
  Tcl_PrintDouble(interp, mmm2d_params.far_cut, buffer);
  Tcl_AppendResult(interp, " ", buffer,(char *) NULL);

  if (mmm2d_params.dielectric_contrast_on) {
    Tcl_PrintDouble(interp, mmm2d_params.delta_mid_top, buffer);
    Tcl_AppendResult(interp, " dielectric-contrasts ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, mmm2d_params.delta_mid_bot, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  }

  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_mmm2d(Tcl_Interp * interp, int argc, char ** argv)
{
  int err;
  double maxPWerror;
  double far_cut = -1;
  double top = 1, mid = 1, bot = 1;
  double delta_top = 0, delta_bot = 0;

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb mmm2d <maximal pairwise error> "
		     "{<fixed far cutoff>} {dielectric <e1> <e2> <e3>} | {dielectric-contrasts <d1> <d2>}", (char *) NULL);
    return TCL_ERROR;
  }
  
  if (! ARG0_IS_D(maxPWerror))
    return TCL_ERROR;
  --argc; ++argv;
  
  if (argc >= 1) {
    if (ARG0_IS_D(far_cut)){
      --argc; ++argv;
    } else {
      Tcl_ResetResult(interp);
    }
  }
  
  if (argc != 0) {
    if (argc == 4 && ARG0_IS_S("dielectric")) {
      if (!ARG_IS_D(1,top) || !ARG_IS_D(2,mid) || !ARG_IS_D(3,bot))
	return TCL_ERROR;

      delta_top = (mid - top)/(mid + top);
      delta_bot = (mid - bot)/(mid + bot);
    }
    else if (argc == 3 && ARG0_IS_S("dielectric-contrasts")) {
      if (!ARG_IS_D(1,delta_top) || !ARG_IS_D(2,delta_bot))
	return TCL_ERROR;
    } else {
      Tcl_AppendResult(interp, "wrong # arguments: inter coulomb mmm2d <maximal pairwise error> "
		       "{<fixed far cutoff>} {dielectric <e1> <e2> <e3>} | {dielectric-contrasts <d1> <d2>}", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if (cell_structure.type != CELL_STRUCTURE_NSQUARE &&
      cell_structure.type != CELL_STRUCTURE_LAYERED) {
    Tcl_AppendResult(interp, "MMM2D requires layered of nsquare cell structure", (char *)NULL);
    return TCL_ERROR;
  }

  if ((err = MMM2D_set_params(maxPWerror, far_cut, delta_top, delta_bot)) > 0) {
    Tcl_AppendResult(interp, mmm2d_errors[err], (char *)NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}


#endif
