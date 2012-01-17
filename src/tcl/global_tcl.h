/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#include "config.h"
#include <tcl.h>

/// Implements the Tcl command setmd. It allows to modify simulation parameters
int tclcommand_setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Implements the Tcl command code_info.  It provides information on the
    Version, Compilation status and the debug status of the used
    code. */
int tclcommand_code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

#endif
