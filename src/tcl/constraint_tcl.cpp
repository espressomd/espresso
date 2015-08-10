/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
/** \file constraint.cpp
    Implementation of \ref constraint.hpp "constraint.h", here it's just the parsing stuff.
*/

#include "config.hpp"

#ifdef CONSTRAINTS

#include <limits>
#include <list>
#include "communication.hpp"
#include "parser.hpp"
#include "TclOutputHelper.hpp"

#include "constraints/ConstraintList.hpp"
#include "constraints/InteractionConstraint.hpp"
#include "shapes/Shape.hpp"
#include "shapes/ShapeList.hpp"
#include "shape_tcl.hpp"

static int tclprint_to_result_Constraint(Tcl_Interp *interp, int i)
{
  return (TCL_OK);
}

int tclcommand_constraint_print(Tcl_Interp *interp)
{
  return (TCL_OK);
}

static void tclprint_to_result_ConstraintForce(Tcl_Interp *interp, int con)
{
  double f[3];
  char buffer[TCL_DOUBLE_SPACE];

  mpi_get_constraint_force(con, f);

  Tcl_PrintDouble(interp, f[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
}


static void tclprint_to_result_n_constraints(Tcl_Interp *interp)
{
  char buffer[TCL_INTEGER_SPACE];
  sprintf(buffer, "%ld", Constraints::list.size());
  Tcl_AppendResult(interp, buffer, (char *) NULL);
}

#endif


static int parse_constraint(Tcl_Interp *interp, int argc, char **argv) {
  std::list<std::string> args;
/* Make a vector from argv and drop 'constraint' from the front */
  args.insert(args.begin(), argv + 1, argv + argc);

  if(args.front() == "mindist_position") {
    args.pop_front();
    return print_mindist_position(interp, args);      
  }
}

int tclcommand_constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef CONSTRAINTS
  int status;
  unsigned int c_num;

  parse_constraint(interp, argc, argv);

  if (argc < 2) return tclcommand_constraint_print(interp);

  if(!strncmp(argv[1], "mindist_position", strlen(argv[1]))) {
    return tclcommand_constraint_mindist_position(interp, argc-2, argv+2);
  }
  if(!strncmp(argv[1], "mindist_position_vec", strlen(argv[1]))) {
    return tclcommand_constraint_mindist_position_vec(interp, argc-2, argv+2);
  }
  else if(!strncmp(argv[1], "wall", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "rhomboid", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "rod", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "plate", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "maze", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "pore", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "slitpore", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "stomatocyte", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "hollow_cone", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "ext_magn_field", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "plane cell", strlen(argv[1]))) {

  }
  else if(!strncmp(argv[1], "force", strlen(argv[1]))) {
    if(argc < 3) {
      Tcl_AppendResult(interp, "which particles force?",(char *) NULL);
      return (TCL_ERROR);
    }
    if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
    if(c_num < 0 || c_num >= Constraints::list.size()) {
      Tcl_AppendResult(interp, "constraint does not exist",(char *) NULL);
      return (TCL_ERROR);
    }
    tclprint_to_result_ConstraintForce(interp, c_num);
    status  = TCL_OK;
  }
  else if(!strncmp(argv[1], "n_constraints", strlen(argv[1]))) {
    tclprint_to_result_n_constraints(interp);
    status=TCL_OK;
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all */
      mpi_bcast_constraint(-2);
      status = TCL_OK;
    }
    else {
      if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
      if(c_num < 0 || c_num >= n_constraints) {
	Tcl_AppendResult(interp, "Can not delete non existing constraint",(char *) NULL);
	return (TCL_ERROR);
      }
      mpi_bcast_constraint(c_num);
      status = TCL_OK;    
    }
  }
  else if (argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    tclprint_to_result_Constraint(interp, c_num);
    status = TCL_OK;
  }
  else {
  //ER "ext_magn_field" was put in the next line //end ER
    Tcl_AppendResult(interp, "possible constraints: wall sphere cylinder maze pore stomatocyte hollow_cone ext_magn_field or constraint delete {c} to delete constraint(s)",(char *) NULL);
    return (TCL_ERROR);
  }

  return gather_runtime_errors(interp, status);

#else /* !defined(CONSTRAINTS) */
  Tcl_AppendResult(interp, "Constraints not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}

