// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file cells.c
 *
 *  This file contains everything related to the link cell
 *  algorithm. 
 *
 *  For more information on cells,
 *  see \ref cells.h "cells.h"
 *   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cells.h"
#include "config.h"
#include "debug.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "communication.h"
#include "utils.h"
#include "verlet.h"
#include "ghosts.h"
#include "domain_decomposition.h"
#include "parser.h"

/* Variables */

/** list of all cells containing particles physically on the local node */
CellPList local_cells = { NULL, 0, 0 };
/** list of all cells containing ghosts */
CellPList ghost_cells = { NULL, 0, 0 };

/** Type of cell structure in use */
CellStructure cell_structure;

/************************************************************/
static void topology_init(int cs, CellPList *local) {
  switch (cs) {
  case CELL_STRUCTURE_CURRENT:
    cell_structure.topology_init(local);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init(local);
    break;
  default:
    fprintf(stderr, "ERROR: attempting to sort the particles in an unknown way\n");
    errexit();
  }
}

/************************************************************/
void cells_init()
{
  /* here local_cells should be empty still */
  dd_topology_init(&local_cells);
}

/************************************************************/
static void free_cellplist(CellPList *cl)
{
  int i;
  for (i = 0; i < cl->n; i++) {
    realloc_particles(&cl->cell[i]->pList, 0);
    free(cl->cell[i]->nList);
    free(cl->cell[i]);
  }
  realloc_cellplist(cl, 0);
}

/************************************************************/
void cells_re_init(int new_cs)
{
  CellPList tmp_local;
  cell_structure.topology_release();
  /* transfer old particle cells to tmp buffer */
  tmp_local = local_cells;
  /* clear local cells */
  init_cellplist(&local_cells);
  /* deallocate the ghosts */
  free_cellplist(&ghost_cells);
  topology_init(new_cs, &tmp_local);
  /* finally deallocate the old cells */
  free_cellplist(&tmp_local);
}

/*************************************************/
int cells_get_n_particles()
{
  int c, cnt = 0;
  for (c = 0; c < local_cells.n; c++)
    cnt += local_cells.cell[c]->pList.n;
  return cnt;
}

/*************************************************/
int cellsystem(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv)
{
  if (argc < 1) {
    Tcl_AppendResult(interp, "usage: cellsystem <system> <params>", (char *)NULL);
    return TCL_ERROR;
  }
  if (ARG0_IS_S("domain_decomposition"))
    mpi_bcast_cell_structure(CELL_STRUCTURE_DOMDEC);
  else {
    Tcl_AppendResult(interp, "unkown cell structure type \"", argv[0],"\"", (char *)NULL);
    return TCL_ERROR;
  }
  return TCL_OK;
}

/*************************************************/
void print_particle_positions()
{
  Cell *cell;
  int c,i,np,cnt=0;
  Particle *part;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->pList.part;
    np   = cell->pList.n;
    for(i=0 ; i < np; i++) {
      fprintf(stderr,"%d: local cell %d contains part id=%d pos=(%f,%f,%f)\n",
	      this_node, c, part[i].p.identity,
	      part[i].r.p[0], part[i].r.p[1], part[i].r.p[2]);
      cnt++;
    }
  }
  fprintf(stderr,"%d: found %d particles\n",this_node,cnt);
}

/*************************************************/
void print_ghost_positions()
{
  Cell *cell;
  int c,i,np,cnt=0;
  Particle *part;

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    part = cell->pList.part;
    np   = cell->pList.n;
    for(i=0 ; i < np; i++) {
      fprintf(stderr,"%d: local cell %d contains ghost id=%d pos=(%f,%f,%f)\n",
	      this_node, c, part[i].p.identity,
	      part[i].r.p[0], part[i].r.p[1], part[i].r.p[2]);
      cnt++;
    }
  }
  fprintf(stderr,"%d: found %d ghosts\n",this_node,cnt);
}
