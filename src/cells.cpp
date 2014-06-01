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
/** \file cells.cpp
 *
 *  This file contains functions for the cell system.
 *
 *  For more information on cells, see cells.hpp
 *   */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "utils.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "integrate.hpp"
#include "initialize.hpp"
#include "communication.hpp"
#include "verlet.hpp"
#include "ghosts.hpp"
#include "domain_decomposition.hpp"
#include "nsquare.hpp"
#include "layered.hpp"

/* Variables */

/** list of all cells. */
Cell *cells = NULL;
/** size of \ref cells */
int n_cells = 0;
/** list of pointers to all cells containing particles physically on the local node. */
CellPList local_cells = { NULL, 0, 0 };
/** list of pointers to all cells containing ghosts. */
CellPList ghost_cells = { NULL, 0, 0 };

/** Type of cell structure in use */
CellStructure cell_structure = { CELL_STRUCTURE_NONEYET };

double max_range = 0.0;

int rebuild_verletlist = 0;

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Extensive Debug function to check the consistency of the cells and
    the particles therein. Use with care! */
void check_cells_consistency()
{
  int c, index;
  IntList used;
  alloc_intlist(&used, n_cells);
  memset(used.e, 0, n_cells*sizeof(int));
  
  for (c = 0; c < local_cells.n; c++) {
    index = (char *)local_cells.cell[c] - (char *)cells;
    if ((index % sizeof(Cell)) != 0) {
      fprintf(stderr, "%d: local cell pointer not even aligned, certainly wrong (local_cell[%d], index=%d).\n", this_node, c, index);
      errexit();
    }
    index /= sizeof(Cell);
    if (index < 0 || index >= n_cells) {
      fprintf(stderr, "%d: local cell pointer out of range, maybe old leftover (local_cell[%d]).\n", this_node, c);
      errexit();
    }
    if (used.e[index]) {
      fprintf(stderr, "%d: local cell is already pointed to (local_cell[%d]).\n", this_node, c);
      errexit();
    }
    used.e[index] = 1;
  }

  for (c = 0; c < ghost_cells.n; c++) {
    index = (char *)ghost_cells.cell[c] - (char *)cells;
    if ((index % sizeof(Cell)) != 0) {
      fprintf(stderr, "%d: ghost cell pointer not even aligned, certainly wrong (ghost_cell[%d], index=%d).\n", this_node, c, index);
      errexit();
    }
    index /= sizeof(Cell);
    if (index < 0 || index >= n_cells) {
      fprintf(stderr, "%d: ghost cell pointer out of range, maybe old leftover (ghost_cell[%d]).\n", this_node, c);
      errexit();
    }
    if (used.e[index]) {
      fprintf(stderr, "%d: ghost cell is already pointed to (ghost_cell[%d]).\n", this_node, c);
      errexit();
    }
    used.e[index] = 1;
  }
  for (c = 0; c < n_cells; c++)
    if (!used.e[c]) {
      fprintf(stderr, "%d: cell %d is not used anywhere.\n", this_node, c);
      errexit();
    }
  realloc_intlist(&used, 0);
}

/** Switch for choosing the topology release function of a certain
    cell system. */
static void topology_release(int cs) {
  switch (cs) {
  case CELL_STRUCTURE_NONEYET:
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_release(cell_structure.type);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_release();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_release();
    break;
  case CELL_STRUCTURE_LAYERED:
    layered_topology_release();
    break;
  default:
    fprintf(stderr, "INTERNAL ERROR: attempting to sort the particles in an unknown way (%d)\n", cs);
    errexit();
  }
}

/** Switch for choosing the topology init function of a certain
    cell system. */
static void topology_init(int cs, CellPList *local) {
  switch (cs) {
  case CELL_STRUCTURE_NONEYET:
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_init(cell_structure.type, local);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init(local);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_init(local);
    break;
  case CELL_STRUCTURE_LAYERED:
    layered_topology_init(local);
    break;
  default:
    fprintf(stderr, "INTERNAL ERROR: attempting to sort the particles in an unknown way (%d)\n", cs);
    errexit();
  }
}

/*@}*/

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

/************************************************************/

void cells_re_init(int new_cs)
{
  CellPList tmp_local;
  Cell *tmp_cells;
  int tmp_n_cells,i;

  CELL_TRACE(fprintf(stderr, "%d: cells_re_init: convert type (%d->%d)\n", this_node, cell_structure.type, new_cs));

  invalidate_ghosts();

  /* 
     CELL_TRACE({
     int p;
     for (p = 0; p < n_part; p++)
     if (local_particles[p])
     fprintf(stderr, "%d: cells_re_init: got particle %d\n", this_node, p);
     }
     );
  */

  topology_release(cell_structure.type);
  /* MOVE old local_cell list to temporary buffer */
  memcpy(&tmp_local,&local_cells,sizeof(CellPList));
  init_cellplist(&local_cells);

  /* MOVE old cells to temporary buffer */
  tmp_cells   = cells;
  tmp_n_cells = n_cells;
  cells   = NULL;
  n_cells = 0;

  topology_init(new_cs, &tmp_local);

  particle_invalidate_part_node();

  /* finally deallocate the old cells */
  realloc_cellplist(&tmp_local,0);
  for(i=0;i<tmp_n_cells;i++)
    realloc_particlelist(&tmp_cells[i],0);

  free(tmp_cells);
  CELL_TRACE(fprintf(stderr, "%d: old cells deallocated\n",this_node));

  /*
    CELL_TRACE({
    int p;
    for (p = 0; p < n_part; p++)
    if (local_particles[p])
    fprintf(stderr, "%d: cells_re_init: now got particle %d\n", this_node, p);
    }
    );
  */

  /* to enforce initialization of the ghost cells */
  resort_particles = 1;

#ifdef ADDITIONAL_CHECKS
  check_cells_consistency();
#endif

  on_cell_structure_change();
}

/************************************************************/

void realloc_cells(int size)
{
  int i;
  CELL_TRACE(fprintf(stderr, "%d: realloc_cells %d\n", this_node, size));
  /* free all memory associated with cells to be deleted. */
  for(i=size; i<n_cells; i++) {
    realloc_particlelist(&cells[i],0);
  }
  /* resize the cell list */
  if(size != n_cells) {
    cells = (Cell *) realloc(cells, sizeof(Cell)*size);
  }
  /* initialize new cells */
  for(i=n_cells; i<size; i++) {
    init_particlelist(&cells[i]);
  }
  n_cells = size;
}  

/*************************************************/

void announce_resort_particles()
{
  int sum;
  
  MPI_Allreduce(&resort_particles, &sum, 1, MPI_INT, MPI_SUM, comm_cart);
  resort_particles = (sum > 0) ? 1 : 0;
  
  INTEG_TRACE(fprintf(stderr,"%d: announce_resort_particles: resort_particles=%d\n",
		      this_node, resort_particles));
}

/*************************************************/

int cells_get_n_particles()
{
  int c, cnt = 0;
  for (c = 0; c < local_cells.n; c++)
    cnt += local_cells.cell[c]->n;
  return cnt;
}

/*************************************************/

void print_local_particle_positions()
{
  Cell *cell;
  int c,i,np,cnt=0;
  Particle *part;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
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

#ifdef CELL_DEBUG

static void dump_particle_ordering()
{
  /* Loop local cells */
  for (int c = 0; c < local_cells.n; c++) {
    Cell *cell  = local_cells.cell[c];
    Particle *p = cell->part;
    int np      = cell->n;

    fprintf(stderr, "%d: cell %d has particles", this_node, c);

    /* Loop cell particles */
    for(int i=0; i < np; i++) {
      fprintf(stderr, " %d", p[i].p.identity);
    }
    fprintf(stderr, "\n");
  }
}

#endif // CELL_TRACE

/*************************************************/

void cells_resort_particles(int global_flag)
{
  CELL_TRACE(fprintf(stderr, "%d: entering cells_resort_particles %d\n", this_node, global_flag));

  invalidate_ghosts();

  particle_invalidate_part_node();
  n_verlet_updates++;

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_exchange_and_sort_particles(global_flag);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_balance_particles(global_flag);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_exchange_and_sort_particles(global_flag);
    break;
  }

#ifdef ADDITIONAL_CHECKS
  /* at the end of the day, everything should be consistent again */
  check_particle_consistency();
#endif

  ghost_communicator(&cell_structure.ghost_cells_comm);
  ghost_communicator(&cell_structure.exchange_ghosts_comm);

  resort_particles = 0;
  rebuild_verletlist = 1;

  on_resort_particles();

  CELL_TRACE(dump_particle_ordering());
  CELL_TRACE(fprintf(stderr, "%d: leaving cells_resort_particles\n", this_node));
}

/*************************************************/

static int compare_particles(const void *a, const void *b)
{ 
  int id_a = static_cast<const Particle *>(a)->p.identity;
  int id_b = static_cast<const Particle *>(b)->p.identity;
  return id_a - id_b;
}

void local_sort_particles()
{
  CELL_TRACE(fprintf(stderr, "%d: entering local_sort_particles\n", this_node));

  /* first distribute strictly on nodes */
  cells_resort_particles(CELL_GLOBAL_EXCHANGE);

  CELL_TRACE(fprintf(stderr, "%d: sorting local cells\n", this_node));

  /* now sort the local cells */
  for (int c = 0; c < local_cells.n; c++) {
    Cell *cell  = local_cells.cell[c];
    Particle *p = cell->part;
    int np      = cell->n;

#ifdef CELL_DEBUG
    for (int id = 0; id < np; ++id) {
      Cell *tgt_cell = cell_structure.position_to_cell(p[id].r.p);
      if (tgt_cell != cell) {
        fprintf(stderr, "%d: particle %d at position %lf %lf %lf is not in its expected cell. Have %ld, expected %ld\n", this_node, p[id].p.identity, p[id].r.p[0], p[id].r.p[1], p[id].r.p[2], (cell - cells)/sizeof(Cell*), (tgt_cell - cells) /sizeof(Cell*)); 
      }
    }
#endif

    qsort(p, np, sizeof(Particle), compare_particles);
    update_local_particles(cell);
  }

  CELL_TRACE(dump_particle_ordering());
  CELL_TRACE(fprintf(stderr, "%d: leaving local_sort_particles\n", this_node));
}

/*************************************************/

void cells_on_geometry_change(int flags)
{
  if (max_cut > 0.0) {
    if (skin >= 0.0)
      max_range = max_cut + skin;
    else
      /* if the skin is not yet set, assume zero. */
      max_range = max_cut;
  }
  else
    /* if no interactions yet, we also don't need a skin */
    max_range = 0.0;

  CELL_TRACE(fprintf(stderr,"%d: on_geometry_change with max range %f\n", this_node, max_range));

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    dd_on_geometry_change(flags);
    break;
  case CELL_STRUCTURE_LAYERED:
    /* there is no fast version, always redo everything. */
    cells_re_init(CELL_STRUCTURE_LAYERED);
    break;
  case CELL_STRUCTURE_NSQUARE:
    /* this cell system doesn't need to react, just tell
       the others */
    on_boxl_change();
    break;
  }
}

/*************************************************/

void check_resort_particles()
{
  int i, c, np;
  Cell *cell;
  Particle *p;
  double skin2 = SQR(skin/2.0);

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      /* Verlet criterion check */
      if(distance2(p[i].r.p, p[i].l.p_old) > skin2) resort_particles = 1;
    }
  }
  announce_resort_particles();
}

/*************************************************/
void cells_update_ghosts()
{
  /* if dd.use_vList is set, it so far means we want EXACT sorting of the particles.*/
  if (dd.use_vList == 0)
    resort_particles = 1;

  if (resort_particles) {
    /* Communication step:  number of ghosts and ghost information */
    cells_resort_particles(CELL_NEIGHBOR_EXCHANGE);
  }
  else
    /* Communication step: ghost information */
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
}

/*************************************************/

void print_ghost_positions()
{
  Cell *cell;
  int c,i,np,cnt=0;
  Particle *part;

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0 ; i < np; i++) {
      fprintf(stderr,"%d: local cell %d contains ghost id=%d pos=(%f,%f,%f)\n",
	      this_node, c, part[i].p.identity,
	      part[i].r.p[0], part[i].r.p[1], part[i].r.p[2]);
      cnt++;
    }
  }
  fprintf(stderr,"%d: found %d ghosts\n",this_node,cnt);
}
