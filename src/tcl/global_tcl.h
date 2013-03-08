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
#ifndef GLOBAL_TCL_H
#define GLOBAL_TCL_H
#include "parser.h"

/// Implements the Tcl command setmd. It allows to modify simulation parameters
int tclcommand_setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/// Function signature for the write callback procedures for global variables
typedef int (SetCallback)(Tcl_Interp *interp, void *data);

/// Register a global variable, should be used only in \ref register_global_variables
void register_global_callback(int field, SetCallback *callback);

#endif
