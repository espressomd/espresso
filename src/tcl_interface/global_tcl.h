/*
  Copyright (C) 2010,2011 The ESPResSo project
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
/** \file global.h
    This file contains the code for access to globally defined
    variables using the script command setmd. \ref add_vars "Here"
    you can find details on how to add new variables in the interpreter's
    space.
*/

#include <tcl.h>

/**********************************************
 * misc procedures
 **********************************************/

/// Implements the Tcl command setmd. It allows to modify simulation parameters
int tclcommand_setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Implements the Tcl command code_info.  It provides information on the
    Version, Compilation status and the debug status of the used
    code. */
int tclcommand_code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

#endif
