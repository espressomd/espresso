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

/* Variables */

/** list of all cells containing particles physically on the local node */
CellPList local_cells = { NULL, 0, 0 };
/** list of all cells containing ghosts */
CellPList ghost_cells = { NULL, 0, 0 };

/** Type of cell structure in use */
CellStructure cell_structure;

/************************************************************/
void cells_init()
{
  dd_cells_init(&local_cells);
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
void cells_re_init() 
{
  CellPList tmp_local;
  cell_structure.topology_release();
  /* transfer old particle cells to tmp buffer */
  tmp_local = local_cells;
  /* clear local cells */
  init_cellplist(&local_cells);
  /* deallocate the ghosts */
  free_cellplist(&ghost_cells);
  cell_structure.topology_init(&tmp_local);
  /* finally deallocate the old cells */
  free_cellplist(&tmp_local);
}

#if 0
  int i,j,ind;
  int old_n_cells, old_ghost_cell_grid[3];
  Cell         *old_cells;
  ParticleList *pl;
  Particle     *part;

#ifdef ADDITIONAL_CHECKS
  int part_cnt_old, part_cnt_new;
#endif

  /* first move particles to their nodes. Necessary if
     box length has changed. May also be called in script mode,
     so better also invalidate the node pointers. */
  particle_invalidate_part_node();
  invalidate_ghosts();
  exchange_and_sort_part();

  CELL_TRACE(fprintf(stderr,"%d: cells_re_init \n",this_node));

  /* 1: store old cell grid */
  old_cells = cells;
  old_n_cells = n_cells;
  for(i=0;i<3;i++) old_ghost_cell_grid[i] = ghost_cell_grid[i];
 
  /* 2: setup new cell grid */
  /* 2a: set up dimensions of the cell grid */
  calc_cell_grid();  
  /* 2b: there should be a reasonable number of cells only!
     But we will deal with that later... */
  /* 2c: allocate new cell structure */
  cells  = (Cell *)malloc(n_cells*sizeof(Cell));
  /* 2d: allocate particle arrays */
  for(i=0; i<n_cells; i++) init_cell(&cells[i]);
  /* 2e: init cell neighbors */
  for(i=0; i<n_cells; i++) init_cell_neighbors(i);
 
  /* 3: Transfer Particle data from old to new cell grid */
  for(i=0;i<old_n_cells;i++) {
    pl = &(old_cells[i].pList);
    if(is_inner_cell(i,old_ghost_cell_grid)) {
      for(j=0; j<pl->n; j++) {
	part = &(pl->part[j]);
	ind  = pos_to_cell_grid_ind(part->r.p);
	append_unindexed_particle(&(cells[ind].pList),part);
      }
    }
    if(pl->max>0) free(pl->part);
    if(old_cells[i].n_neighbors>0) {
      for(j=0; j<old_cells[i].n_neighbors; j++) free(old_cells[i].nList[j].vList.pair);
      free(old_cells[i].nList);
    }
  }

  for(i=0;i<n_cells;i++) {
    if(is_inner_cell(i,ghost_cell_grid))
      update_local_particles(&(cells[i].pList));
  }

#ifdef ADDITIONAL_CHECKS
  /* check particle transfer */
  part_cnt_old=0;
  for(i=0;i<old_n_cells;i++) 
    if(is_inner_cell(i,old_ghost_cell_grid)) 
      part_cnt_old += old_cells[i].pList.n;
  part_cnt_new=0;
  for(i=0;i<n_cells;i++) 
    if(is_inner_cell(i,ghost_cell_grid)) 
      part_cnt_new += cells[i].pList.n;
  if(part_cnt_old != part_cnt_new) 
    CELL_TRACE(fprintf(stderr,"%d: cells_re_init: lost particles: old grid had %d new grid has %d particles.\n",this_node,part_cnt_old, part_cnt_new));
#endif

  CELL_TRACE(fprintf(stderr,"%d: cell_grid (%d %d %d) \n",this_node,cell_grid[0],cell_grid[1],cell_grid[2]));

  free(old_cells);
  /* cell structure initialized. */
  rebuild_verletlist = 1;
}

#endif

/*************************************************/
void cells_changed_topology()
{
  cells_re_init();
}

/*************************************************/
int cells_get_n_particles()
{
  int c, cnt = 0;
  for (c = 0; c < local_cells.n; c++)
    cnt += local_cells.cell[c]->pList.n;
  return cnt;
}

#if 0
/*************************************************/
Particle *cells_alloc_particle(int id, double pos[3])
{
  int rl;
  int ind = pos_to_cell_grid_ind(pos);
  ParticleList *pl = &cells[ind].pList;
  Particle *pt;

  pl->n++;
  rl = realloc_particles(pl, pl->n);
  pt = &pl->part[pl->n - 1];
  init_particle(pt);

  pt->p.identity = id;
  memcpy(pt->r.p, pos, 3*sizeof(double));
  if (rl)
    update_local_particles(&cells[ind].pList);
  else
    local_particles[pt->p.identity] = pt;

  return pt;
}
#endif

/************************************************************/
/*******************  privat functions  *********************/
/************************************************************/

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
