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
/** \file config_tcl.c
 *
 *  contains code_info and version stuff.
*/
#include "utils.h"
#include <tcl.h>

static int tclcommand_code_info_version(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, PACKAGE_NAME, "-", ESPRESSO_VERSION, (char *) NULL);
  return (TCL_OK);
}

/** callback for compilation status. */
static int tclcommand_code_info_compilation(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Compilation status ", (char *) NULL);
  for (int i=0; i < NUM_FEATURES; i++) {
    Tcl_AppendResult(interp, "{ ", (char *) NULL);
    Tcl_AppendResult(interp, FEATURES[i], (char *) NULL);
    Tcl_AppendResult(interp, " } ", (char *) NULL);
  }
  Tcl_AppendResult(interp, "}", (char *) NULL);
  return (TCL_OK);
}

int tclcommand_code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  if (argc < 2) {
    tclcommand_code_info_version(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    tclcommand_code_info_compilation(interp);
  }
  else {
    if(!strncmp(argv[1], "version" , strlen(argv[1]))) {
      tclcommand_code_info_version(interp);
    }
    else if(!strncmp(argv[1], "compilation" , strlen(argv[1]))) {
      tclcommand_code_info_compilation(interp);
    }
    else if(!strncmp(argv[1], "debug" , strlen(argv[1]))) {
      tclcommand_code_info_compilation(interp);
    }
    else {
      Tcl_AppendResult(interp, "info ",argv[1]," not known!", (char *) NULL);
    }
  }
  return (TCL_OK);
}
