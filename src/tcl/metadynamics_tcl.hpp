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

#ifndef _METADYNAMICS_TCL_H
#define _METADYNAMICS_TCL_H
#include "parser.hpp"

#ifdef METADYNAMICS

/*********************************
 * functions
 *********************************/
/** Implementation of the Tcl command \ref tclcommand_metadynamics. This function
 *  allows to change the parameters of metadynamics */
int tclcommand_metadynamics(ClientData data, Tcl_Interp *interp, int argc, char **argv);
/** Print metadynamics options and parameters */
int tclcommand_metadynamics_print_status(Tcl_Interp *interp);
int tclcommand_metadynamics_print_usage(Tcl_Interp *interp, int argc, char **argv);
int tclcommand_metadynamics_parse_off(Tcl_Interp *interp, int argc, char **argv);
/** Reaction coordinates available */
int tclcommand_metadynamics_parse_distance(Tcl_Interp *interp, int argc, char **argv);
int tclcommand_metadynamics_parse_relative_z(Tcl_Interp *interp, int argc, char **argv);
/** Input/Output stuff */
int tclcommand_metadynamics_print_stat(Tcl_Interp *interp, int argc, char **argv);
int tclcommand_metadynamics_parse_load_stat(Tcl_Interp *interp, int argc, char **argv);

#endif

#endif
