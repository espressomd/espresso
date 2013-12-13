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
#ifndef _LBBOUNDARIES_TCL_H
#define _LBBOUNDARIES_TCL_H

#include "parser.hpp"
#include "lb-boundaries.hpp"


#if defined (LB_BOUNDARIES) || defined (LB_BOUNDARIES_GPU)

// TCL Parser functions
LB_Boundary *generate_lbboundary();
int tclcommand_lbboundary(ClientData _data, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_wall(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_sphere(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_cylinder(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_rhomboid(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_pore(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_stomatocyte(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_hollow_cone(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_lbboundary_box(LB_Boundary *lbb, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_printLbBoundaryToResult(Tcl_Interp *interp, int i);

#endif /* LB_BOUNDARIES || LB_BOUNDARIES_GPU */

#endif /* _LBBOUNDARIES_TCL_H */
