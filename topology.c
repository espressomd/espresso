// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.

/** \file topology.c
 *
 *  This file contains functions for handling the system topology.
 *
 *  For more information see \ref topology.h "topology.h"
 *   */

#include "utils.h"
#include "parser.h"
#include "topology.h"
#include "statistics_chain.h"

int     n_molecules = -1;
Molecule *molecules = NULL;

void realloc_molecules(int size)
{
  int m;
  for(m = size; m < n_molecules; m++)
    realloc_intlist(&molecules[m].part, 0);

  molecules = realloc(molecules, size*sizeof(Molecule));

  if (n_molecules < 0)
    n_molecules = 0;
  for(m = n_molecules; m < size; m++)
    init_intlist(&molecules[m].part);

  n_molecules = size;
}

int print_structure_info(Tcl_Interp *interp)
{
  char buffer[TCL_INTEGER_SPACE + 1];
  int m, i;
  for (m = 0; m < n_molecules; m++) {
    sprintf(buffer, "%d ", molecules[m].type);
    Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
    for (i = 0; i < molecules[m].part.n; i++) {
      sprintf(buffer, "%d ", molecules[m].part.e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);      
    }
    Tcl_AppendResult(interp, "} ", (char *)NULL);
  }
  return TCL_OK;
}

int parse_generic_structure_info(Tcl_Interp *interp, int argc, char **argv)
{
  int arg;
  IntList il;

  init_intlist(&il);

  realloc_molecules(argc);

  for (arg = 0; arg < argc; arg++) {
    if (!ARG_IS_INTLIST(arg, il)) {
      realloc_molecules(0);
      realloc_intlist(&il, 0);
      return TCL_ERROR;
    }
    molecules[arg].type = il.e[0];
    realloc_intlist(&molecules[arg].part, molecules[arg].part.n = il.n - 1);
    memcpy(molecules[arg].part.e, &il.e[1], (il.n - 1)*sizeof(int));
  }
  realloc_intlist(&il, 0);
  return TCL_OK;
}

int parse_analyze_set_topology(Tcl_Interp *interp, int argc, char **argv)
{
  if (argc == 0)
    return print_structure_info(interp);

  if (ARG0_IS_S("chains"))
    return parse_chain_structure_info(interp, argc - 1, argv + 1);

  return parse_generic_structure_info(interp, argc, argv);
}

