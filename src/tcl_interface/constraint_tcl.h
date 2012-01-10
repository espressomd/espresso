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
#ifndef CONSTRAINT_TCL_H
#define CONSTRAINT_TCL_H
#include "statistics.h"
#include "energy.h"
#include "forces.h"
#include "grid.h"
#include "errorhandling.h"
#include "tunable_slip.h"
#include <tcl.h>

/** \file constraint_tcl.h
 */

#ifdef CONSTRAINTS

int tclprint_to_result_Constraint(Tcl_Interp *interp, int i);
int tclcommand_constraint_print(Tcl_Interp *interp);
void tclprint_to_result_ConstraintForce(Tcl_Interp *interp, int con);
int tclcommand_constraint_parse_wall(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_sphere(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_cylinder(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_pore(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_rod(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_plate(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_maze(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_ext_magn_field(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_parse_plane_cell(Constraint *con, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint_mindist_position(Tcl_Interp *interp, int argc, char **argv);
int tclcommand_constraint(ClientData _data, Tcl_Interp *interp, int argc, char **argv);

#endif

#endif
