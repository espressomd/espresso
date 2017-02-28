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
/** \file elc_tcl.cpp
 *
 *  Implementation of \ref elc_tcl.hpp
 */
#include "elc_tcl.hpp"
#include "elc.hpp"

#ifdef P3M

int tclprint_to_result_ELC(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, elc_params.maxPWerror, buffer);
  Tcl_AppendResult(interp, "} {coulomb elc ", buffer, (char *) NULL);
  Tcl_PrintDouble(interp, elc_params.gap_size, buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  Tcl_PrintDouble(interp, elc_params.far_cut, buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  if (!elc_params.neutralize)
    Tcl_AppendResult(interp, " noneutralization", (char *) NULL);
  if (elc_params.const_pot_on) {
    Tcl_PrintDouble(interp, elc_params.pot_diff, buffer);
    Tcl_AppendResult(interp, " capacitor ", buffer, (char *) NULL);
  }
  else if (elc_params.dielectric_contrast_on) {
    Tcl_PrintDouble(interp, elc_params.di_mid_top, buffer);
    Tcl_AppendResult(interp, " dielectric-contrasts ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, elc_params.di_mid_bot, buffer);
    Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  }
  return TCL_OK;
}

int tclcommand_inter_coulomb_parse_elc_params(Tcl_Interp * interp, int argc, char ** argv)
{
  double pwerror;
  double gap_size;
  double far_cut = -1;
  double top = 1, mid = 1, bot = 1;
  double delta_top = 0, delta_bot = 0;
  int neutralize = 1;
  double pot_diff = 0;
  int const_pot_on = 0;

  if (argc < 2) {
    Tcl_AppendResult(interp, "either nothing or elc <pwerror> <minimal layer distance> {<cutoff>} <{dielectric <di_top> <di_mid> <di_bottom>} | {dielectric-contrasts <d1> <d2>} | {capacitor <dU>}> {noneutralization} expected, not \"",
		     argv[0], "\"", (char *)NULL);
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

    while (argc > 0) {
      if (ARG0_IS_S("noneutralization") || ARG0_IS_S("-noneutralization")) {
	neutralize = 0;
	argc--; argv++;
      }
      else if (argc >= 4 && ARG0_IS_S("dielectric")) {
  	Tcl_AppendResult(interp, "There seems to be an error when using ELC with dielectric constrasts. If you are sure you want to use it, you have to deactivate this message manually. ", (char *)NULL);
	return TCL_ERROR;
	// just a dummy, not used, as it is only printed for information
	// purposes. We need to calculate it
	double space_layer_dummy;

	if (!ARG_IS_D(1,top) || !ARG_IS_D(2,mid) || !ARG_IS_D(3,bot))
	  return TCL_ERROR;
        delta_top = (mid - top)/(mid + top);
        delta_bot = (mid - bot)/(mid + bot);
	argc -= 4; argv += 4;

	if (argc > 0 && ARG_IS_D(4, space_layer_dummy)) {
	  argc--; argv++;
	}
      }
      else if (argc >= 3 && ARG0_IS_S("dielectric-contrasts")) {
  	Tcl_AppendResult(interp, "There seems to be an error when using ELC with dielectric constrasts. If you are sure you want to use it, you have to deactivate this message manually. ", (char *)NULL);
	return TCL_ERROR;

        if (!ARG_IS_D(1,delta_top) || !ARG_IS_D(2,delta_bot))
	  return TCL_ERROR;
        argc -= 3; argv += 3;
      }
      else if (argc >= 1 && ARG0_IS_S("capacitor")) {
  	Tcl_AppendResult(interp, "There seems to be an error when using ELC with dielectric constrasts. If you are sure you want to use it, you have to deactivate this message manually. ", (char *)NULL);
	return TCL_ERROR;
	if (!ARG_IS_D(1,pot_diff))
	   return TCL_ERROR;
 	argc -= 2; argv += 2;
        const_pot_on = 1;
	delta_top = -1; delta_bot = -1;
      }
      else {
	Tcl_AppendResult(interp, "either nothing or elc <pwerror> <minimal layer distance> {<cutoff>}  <{dielectric <di_top> <di_mid> <di_bottom>} | {dielectric-contrasts <d1> <d2>} | {capacitor <dU>}> {noneutralization} expected, not \"", argv[0], "\"", (char *)NULL);
	return TCL_ERROR;
      }
    }
  }
  CHECK_VALUE(ELC_set_params(pwerror, gap_size, far_cut, neutralize, delta_top, delta_bot, const_pot_on, pot_diff),
	      "choose a 3d electrostatics method prior to ELC");
}

#endif
