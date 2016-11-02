/*
  Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project
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
/** \file lees_edwards_tcl.cpp TCL interface for setting the Lees Edwards PBC shear rate.
 *
 *  For more information on LE-PBCs
 *  see \ref lees_edwards.hpp "less_edwards.hpp". 
*/
#include "utils.hpp"
#include "parser.hpp"
#include "communication.hpp"
#include "global.hpp"
#include "lees_edwards.hpp"

int tclcommand_lees_edwards_offset(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  /* check usage */
#ifndef LEES_EDWARDS
  Tcl_AppendResult(interp, "ERROR: Requested a feature which was not compiled in, please reconfigure and recompile with LEES_EDWARDS defined.", (char *)NULL); 
  return (TCL_ERROR);
#else

  char buffer[120 + TCL_DOUBLE_SPACE];
  double new_offset;


  if (argc < 1 || argc > 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: lees_edwards_offset { new_offset }", (char *)NULL); return (TCL_ERROR);
  }

  /* return the old value, whether we have set a new one or not */
  sprintf(buffer, "%f ", lees_edwards_offset);

  /* set the new value if requested */
  if (argc == 2 ){
    if (Tcl_GetDouble(interp, argv[1], &new_offset) == TCL_ERROR) return (TCL_ERROR);
    lees_edwards_offset = new_offset;
    mpi_bcast_parameter(FIELD_LEES_EDWARDS_OFFSET);
  }
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  return gather_runtime_errors(interp, TCL_OK);
#endif
}



