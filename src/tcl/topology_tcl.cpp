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

/** \file topology.c
 *
 *  This file contains functions for handling the system topology.
 *
 *  For more information see topology.h
 *   */

//#include "utils.hpp"
#include "parser.hpp"
#include "topology.hpp"
#include "statistics_chain_tcl.hpp"
#include "particle_data.hpp"
//#include "cells.hpp"
#include "communication.hpp"
//#include "molforces.hpp"

int tclcommand_analyze_set_print_all(Tcl_Interp *interp)
{
  char buffer[TCL_INTEGER_SPACE + 2];
  int m, i;
  for (m = 0; m < n_molecules; m++) {
    sprintf(buffer, "%d ", topology[m].type);
    Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
    for (i = 0; i < topology[m].part.n; i++) {
      sprintf(buffer, "%d ", topology[m].part.e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);      
    }
    Tcl_AppendResult(interp, "} ", (char *)NULL);
  }
  return TCL_OK;
}

int tclcommand_analyze_parse_generic_structure(Tcl_Interp *interp, int argc, char **argv)
{
  int arg;
  IntList il;
  init_intlist(&il);

  realloc_topology(argc);
  
  for (arg = 0; arg < argc; arg++) {
    if (!ARG_IS_INTLIST(arg, il)) {
      realloc_topology(0);
      realloc_intlist(&il, 0);
      return TCL_ERROR;
    }
    topology[arg].type = il.e[0];
    realloc_intlist(&topology[arg].part, topology[arg].part.n = il.n - 1);
    memcpy(topology[arg].part.e, &il.e[1], (il.n - 1)*sizeof(int));
  }
  realloc_intlist(&il, 0);

  return TCL_OK;
}

int tclcommand_analyze_set_parse_topo_part_sync(Tcl_Interp *interp) {
  int i,j,ntopoparts;

  if (n_molecules <= 0) {
    Tcl_AppendResult(interp, "Can't sync molecules to particle info: No molecules defined ", (char *)NULL);
    return TCL_ERROR;
  }
  if (n_total_particles <= 0) {
    Tcl_AppendResult(interp, "Can't sync molecules to particle info: No particles defined ", (char *)NULL);
    return TCL_ERROR;
  }

  /* Check to see that the number of particles in the topology info
     does not exceed the total number of particles */
  ntopoparts = 0;
  for ( i = 0 ; i < n_molecules ; i ++ ) {
    for ( j = 0 ; j < topology[i].part.n ; j++ ) {
      ntopoparts += 1;
    }
  }
  if ( ntopoparts > n_total_particles ) {
    Tcl_AppendResult(interp, "Can't sync molecules to particle info: Topology contains more particles than actually exist ", (char *)NULL);
    return TCL_ERROR;
  }

  if ( !mpi_sync_topo_part_info()) {
    Tcl_AppendResult(interp, "Error syncronising molecules to particle info", (char *)NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}

int tclcommand_analyze_set_parse_trapmol(Tcl_Interp *interp, int argc, char **argv)
{

#ifdef MOLFORCES
#ifdef EXTERNAL_FORCES
  int trap_flag = 0;
  int noforce_flag =0;
  int i;
#endif
#endif
  int mol_num;
  double spring_constant;
  double drag_constant;
  int isrelative;
  DoubleList trap_center;
  IntList trap_coords;
  IntList noforce_coords;
  char usage[] = "trapmol usage: <mol_id> { <xpos> <ypos> <zpos> } <isrelative> <spring_constant> <drag_constant> coords   { <trapped_coord> <trapped_coord> <trapped_coord> } noforce_coords {<noforce_coord> <noforce_coord> <noforce_coord>}";

  init_doublelist(&trap_center);
  init_intlist(&trap_coords);
  alloc_intlist(&trap_coords,3);
  init_intlist(&noforce_coords);
  alloc_intlist(&noforce_coords,3);
  /* Unless coords are specified the default is just to trap it completely */
  trap_coords.e[0] = 1;
  trap_coords.e[1] = 1;
  trap_coords.e[2] = 1;

  Tcl_ResetResult(interp);
  /* The first argument should be a molecule number */
  if (!ARG0_IS_I(mol_num)) {
    Tcl_AppendResult(interp, "first argument should be a molecule id", (char *)NULL);
    Tcl_AppendResult(interp, usage, (char *)NULL); 
    return TCL_ERROR;
  } else {
    /* Sanity checks */
    if (mol_num > n_molecules) {
      Tcl_AppendResult(interp, "trapmol: cannot trap mol %d because it does not exist",mol_num , (char *)NULL);
    return TCL_ERROR;
    }
    argc--;
    argv++;
  }

  /* The next argument should be a double list specifying the trap center */
  if (!ARG0_IS_DOUBLELIST(trap_center)) {
    Tcl_AppendResult(interp, "second argument should be a double list", (char *)NULL);
    Tcl_AppendResult(interp, usage , (char *)NULL);
    return TCL_ERROR;
  } else {
    argc -= 1;
    argv += 1;
  }

  /* The next argument should be an integer specifying whether the trap is relative (fraction of box_l) or absolute */
  if (!ARG0_IS_I(isrelative)) {
    Tcl_AppendResult(interp, "third argument should be an integer", (char *)NULL);
    Tcl_AppendResult(interp, usage, (char *)NULL);
    return TCL_ERROR;
  } else {
    argc -= 1;
    argv += 1;
  }

  /* The next argument should be the spring constant for the trap */
  if (!ARG0_IS_D(spring_constant)) {
    Tcl_AppendResult(interp, "fourth argument should be a double", (char *)NULL);
    Tcl_AppendResult(interp, usage, (char *)NULL);
    return TCL_ERROR;
  } else {
    argc -= 1;
    argv += 1;
  }

  /* The next argument should be the drag constant for the trap */
  if (!ARG0_IS_D(drag_constant)) {
    Tcl_AppendResult(interp, "fifth argument should be a double", (char *)NULL);
    Tcl_AppendResult(interp, usage, (char *)NULL);
    return TCL_ERROR;
  } else {
    argc -= 1;
    argv += 1;
  }

  /* Process optional arguments */
  while ( argc > 0 ) {    
    if ( ARG0_IS_S("coords") ) {
      if ( !ARG_IS_INTLIST(1,trap_coords) ) {
	Tcl_AppendResult(interp, "an intlist is required to specify coords", (char *)NULL);
	Tcl_AppendResult(interp, usage, (char *)NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else if ( ARG0_IS_S("noforce_coords")) {
      if ( !ARG_IS_INTLIST(1,noforce_coords) ) {
	Tcl_AppendResult(interp, "an intlist is required to specify coords", (char *)NULL);
	Tcl_AppendResult(interp, usage, (char *)NULL);
	return TCL_ERROR;
      }
      argc -= 2;
      argv += 2;
    } else {
      Tcl_AppendResult(interp, "an option is not recognised", (char *)NULL);
      Tcl_AppendResult(interp, usage, (char *)NULL);
      return TCL_ERROR;
    }      
  }

#ifdef MOLFORCES 
#ifdef EXTERNAL_FORCES 
  for (i = 0; i < 3; i++) {
    if (trap_coords.e[i])
      trap_flag |= COORD_FIXED(i);
  
    if (noforce_coords.e[i])
      noforce_flag |= COORD_FIXED(i);
  }
  if (set_molecule_trap(mol_num, trap_flag,&trap_center,spring_constant, drag_constant, noforce_flag, isrelative) == TCL_ERROR) {
    Tcl_AppendResult(interp, "set topology first", (char *)NULL);
    return TCL_ERROR;
  }
#else
    Tcl_AppendResult(interp, "Error: EXTERNAL_FORCES not defined ", (char *)NULL);
    return TCL_ERROR;
#endif
#endif

  realloc_doublelist(&trap_center,0);
  realloc_intlist(&trap_coords,0);
  realloc_intlist(&noforce_coords,0);
  return TCL_OK;
  
}

int tclcommand_analyze_parse_set(Tcl_Interp *interp, int argc, char **argv)
{
  if (argc == 0)
    return tclcommand_analyze_set_print_all(interp);

  if (ARG0_IS_S("chains")) {
    return tclcommand_analyze_set_parse_chain_topology(interp, argc - 1, argv + 1);
  } else if (ARG0_IS_S("topo_part_sync")) {

    return tclcommand_analyze_set_parse_topo_part_sync(interp);
  } else if (ARG0_IS_S("trapmol")) {
#ifndef MOLFORCES
    Tcl_AppendResult(interp, "Attempt to trap molecule, but MOLFORCES was not defined.  Turn on MOLFORCES in config.h", (char *)NULL);
    return TCL_ERROR;
#else
    return tclcommand_analyze_set_parse_trapmol(interp, argc - 1, argv + 1);
#endif
  } 
  return tclcommand_analyze_parse_generic_structure(interp, argc, argv);
}



