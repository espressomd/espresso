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
#ifndef IA_DATA_TCL_H
#define IA_DATA_TCL_H
#include "parser.h"
#include "forcecap_tcl.h"

/************************************************
 * exported functions
 ************************************************/

/** Implementation of the tcl command "inter". This function
    allows the interaction parameters to be modified.
 */
int tclcommand_inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Implementation of the Tcl function constraint. This function
    allows to set and delete constraints.
 */
int tclcommand_constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv);

/** datafield callback for \ref min_global_cut. Sets the minimal cell size. */
int tclcallback_min_global_cut(Tcl_Interp *interp, void *_data);

#endif
