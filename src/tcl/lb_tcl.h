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
#ifndef _LB_TCL_H
#define _LB_TCL_H
#include "parser.h"

/** Parser for the TCL command lbfluid. */
int tclcommand_lbfluid(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Parser for the lbnode command. */
int tclcommand_lbnode(ClientData data, Tcl_Interp *interp, int argc, char **argv);

int tclcommand_lbfluid_print_interpolated_velocity(Tcl_Interp *interp, int argc, char **argv);

int tclcommand_lbnode_extforce_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** lb boundary command. From \ref lb-boundaries_tcl.c */
extern int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp, int argc, char **argv);

#endif /* LB_TCL_H */
