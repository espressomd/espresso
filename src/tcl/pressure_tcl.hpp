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
/** \file pressure.h
    Pressure calculation. Really similar to \ref energy.h "energy.h".
*/

#ifndef _PRESSURE_TCL_H
#define _PRESSURE_TCL_H
#include "parser.hpp"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Callback for setting \ref nptiso_struct::piston */
int tclcallback_npt_piston(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref nptiso_struct::p_ext */
int tclcallback_p_ext(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref nptiso_struct::p_diff */
int tclcallback_npt_p_diff(Tcl_Interp *interp, void *_data);

/** implementation of 'analyze pressure'
    @param interp Tcl interpreter
    @param argc   arguments
    @param argv   arguments
    @param v_comp flag which enables (1) compensation of the velocities required
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT. */
int tclcommand_analyze_parse_and_print_pressure(Tcl_Interp *interp, int v_comp, int argc, char **argv);

/** Implementation of 'analyze bins' */
int tclcommand_analyze_parse_bins(Tcl_Interp *interp, int argc, char **argv);

/** implementation of 'analyze p_IK1' */
int tclcommand_analyze_parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv);

/** implementation of 'analyze stress_tensor' */
int tclcommand_analyze_parse_and_print_stress_tensor(Tcl_Interp *interp, int v_comp, int argc, char **argv);

/** implementation of 'analyse local_stress_tensor */
int tclcommand_analyze_parse_local_stress_tensor(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
