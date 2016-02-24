/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file blockfile_tcl.cpp
    Implements the blockfile command for writing Tcl-formatted data files.
*/

#include <iostream>
#include <vector>

#include "readpdb.hpp"
#include "parser.hpp"

static void usage(Tcl_Interp *interp) {
  Tcl_AppendResult(interp, "usage:", (char *)NULL);
}

int tclcommand_readpdb(ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  char *pdb_file = NULL;
  char *itp_file = NULL;
  int first_id = -1;
  int first_type = 0;
  int type = -1;
  bool fit = false;
  bool lj_internal = false;
  double lj_rel_cutoff = 2.5;
  bool lj_diagonal = false;

  std::vector<PdbLJInteraction> ljinteractions;

  argc--;
  argv++;

  while(argc > 0) {
    if(ARG0_IS_S("pdb_file")) {
      argc--;
      argv++;
      pdb_file = argv[0];
    } else if (ARG0_IS_S("itp_file")) {
      argc--;
      argv++;
      itp_file = argv[0];
    } else if (ARG0_IS_S("type")) {
      argc--;
      argv++;
      if(!ARG0_IS_I(type)) {
	Tcl_AppendResult(interp, "type takes exactly one integer argument.\n", (char *)NULL);
	return TCL_ERROR;
      }
    } else if (ARG0_IS_S("first_id")) {
      argc--;
      argv++;
      if(!ARG0_IS_I(first_id)) {
	Tcl_AppendResult(interp, "first_id takes exactly one integer argument.\n", (char *)NULL);
	return TCL_ERROR;
      }      
    } else if (ARG0_IS_S("first_type")) {
      argc--;
      argv++;
      if(!ARG0_IS_I(first_type)) {
	Tcl_AppendResult(interp, "first_type takes exactly one integer argument.\n", (char *)NULL);
	return TCL_ERROR;
      }            
    } else if (ARG0_IS_S("lj_rel_cutoff")) {
      argc--;
      argv++;
      if(!ARG0_IS_D(lj_rel_cutoff)) {
	return TCL_ERROR;
      }
    } else if (ARG0_IS_S("rescale_box")) {
      fit = true;
    } else if (ARG0_IS_S("lj_internal")) {
      lj_internal = true;
      lj_diagonal = false;
    } else if (ARG0_IS_S("lj_diagonal")) {
      lj_internal = true;
      lj_diagonal = true;
    } else if (ARG0_IS_S("lj_with")) {
      argc--;
      argv++;
      if(argc < 3)
	return TCL_ERROR;
      struct PdbLJInteraction ljia;
      if(!ARG0_IS_I(ljia.other_type)) {
	return TCL_ERROR;
      }                  
      argc--;
      argv++;
      if(!ARG0_IS_D(ljia.epsilon)) {
	return TCL_ERROR;
      }                  
      argc--;
      argv++;
      if(!ARG0_IS_D(ljia.sigma)) {
	return TCL_ERROR;
      }
      ljinteractions.push_back(ljia);
    }
    else {
      usage(interp);
      return TCL_ERROR;
    }
    argc--;
    argv++;
  }  
  if((type < 0) || (first_id < 0) || (pdb_file == NULL)) {
    usage(interp);
    return TCL_ERROR;
  }
  const int n_part = pdb_add_particles_from_file(pdb_file, first_id, type, ljinteractions, lj_rel_cutoff, itp_file, first_type, fit, lj_internal, lj_diagonal);
  if(!n_part) {
    Tcl_AppendResult(interp, "Could not parse pdb file.", (char *)NULL);
    return TCL_ERROR;
  }
  char buffer[32];
  snprintf(buffer, sizeof(buffer), "%d", n_part);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

