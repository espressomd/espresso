// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.

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
#include "particle_data.h"
#include "cells.h"
#include "communication.h"

int     n_molecules = -1;
Molecule *topology = NULL;

void realloc_topology(int size)
{
  int m;

  for(m = size ; m < n_molecules; m++)
    realloc_intlist(&topology[m].part, 0);


  topology = realloc(topology, size*sizeof(Molecule));

  if (n_molecules < 0)
    n_molecules = 0;
  for(m = n_molecules; m < size; m++)
    init_intlist(&topology[m].part);
  
  n_molecules = size;

}

int print_structure_info(Tcl_Interp *interp)
{
  char buffer[TCL_INTEGER_SPACE + 1];
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

int parse_generic_structure_info(Tcl_Interp *interp, int argc, char **argv)
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

// Parallel function for synchronising topology and particle data
void sync_topo_part_info() {
  int i,j;
  Particle* p;
  for ( i = 0 ; i < n_molecules ; i ++ ) {
    for ( j = 0 ; j < topology[i].part.n ; j++ ) {
      p = local_particles[topology[i].part.e[j]];
      if(!p) { /* Do nothing */ } 
      else {
	p->p.mol_id = i;
      }
    }
  }
}


int parse_sync_topo_part_info(Tcl_Interp *interp) {
  
  if (n_molecules <= 0) {
    Tcl_AppendResult(interp, "Can't sync molecules to particle info: No molecules defined ", (char *)NULL);
    return TCL_ERROR;
  }
  if (n_total_particles <= 0) {
    Tcl_AppendResult(interp, "Can't sync molecules to particle info: No particles defined ", (char *)NULL);
    return TCL_ERROR;
  }

  if ( !mpi_sync_topo_part_info()) {
    Tcl_AppendResult(interp, "Error syncronising molecules to particle info", (char *)NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}

int parse_analyze_set_topology(Tcl_Interp *interp, int argc, char **argv)
{
  if (argc == 0)
    return print_structure_info(interp);

  if (ARG0_IS_S("chains")) {
    return parse_chain_structure_info(interp, argc - 1, argv + 1);
  } else if (ARG0_IS_S("topo_part_sync")) {

    return parse_sync_topo_part_info(interp);
  }

  return parse_generic_structure_info(interp, argc, argv);
}



