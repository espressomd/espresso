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
/** \file mdlc_correction_tcl.cpp
 *
 *  Implementation of \ref mdlc_correction_tcl.hpp
 */
#include "mdlc_correction_tcl.hpp"
#include "mdlc_correction.hpp"

#ifdef DIPOLES

int tclprint_to_result_MDLC(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  
  Tcl_PrintDouble(interp, dlc_params.maxPWerror, buffer);
 Tcl_AppendResult(interp, "} {magnetic mdlc ", buffer, (char *) NULL);
  Tcl_PrintDouble(interp, dlc_params.gap_size, buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  Tcl_PrintDouble(interp, dlc_params.far_cut, buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  return TCL_OK;
}
/* ***************************************************************** */

int tclcommand_inter_magnetic_parse_mdlc_params(Tcl_Interp * interp, int argc, char ** argv)
{
  double pwerror;
  double gap_size;
  double far_cut = -1;
 
  MDLC_TRACE(fprintf(stderr, "%d: tclcommand_inter_magnetic_parse_mdlc_params().\n", this_node));
  
  if (argc < 2) {
    Tcl_AppendResult(interp, "either nothing or mdlc <pwerror> <minimal layer distance> {<cutoff>}  expected, not \"", argv[0], "\"", (char *)NULL);
    return TCL_ERROR;
  }
  if (!ARG0_IS_D(pwerror))
    return TCL_ERROR;
  if (!ARG1_IS_D(gap_size))
    return TCL_ERROR;

  argc -= 2; argv += 2;

  if (argc > 0) {
    // if there, parse away manual cutoff
    if(ARG0_IS_D(far_cut)) {
      argc--; argv++;
    }
    else
      Tcl_ResetResult(interp);

    if(argc > 0) {
	Tcl_AppendResult(interp, "either nothing or mdlc <pwerror> <minimal layer distance=size of the gap without particles> {<cutoff>}   expected, not \"", argv[0], "\"", (char *)NULL);
	return TCL_ERROR;
    }
  }

  CHECK_VALUE(mdlc_set_params(pwerror,gap_size,far_cut),
	      "choose a 3d electrostatics method prior to use mdlc");
}
/* ***************************************************************** */

#endif
