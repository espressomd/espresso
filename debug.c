/** \file debug.c
    Implements the malloc replacements as described in \ref debug.h "debug.h". */

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "communication.h"
#include "debug.h"
#include "cells.h"
#include "grid.h"

#if defined FORCE_CORE || defined MPI_CORE
int regular_exit = 0;
#else
int regular_exit = 1;
#endif
static int core_done = 0;

int check_id = 91;

void core()
{
  if (!core_done && !regular_exit) {
    core_done = 1;
    fprintf(stderr, "%d: forcing core dump on irregular exit\n", this_node);
    kill(getpid(), SIGSEGV);
  }
}

void check_particle_consistency()
{
  ParticleList *pl;
  Particle *part;
  Cell *cell;
  int n, np, dir, c, nc;
  int cell_part_cnt=0, ghost_part_cnt=0, local_part_cnt=0;
  int cell_err_cnt=0;

  /* checks: part_id, part_pos, local_particles id */
  for(c=0; c<n_cells; c++) {
    if(is_inner_cell(c,ghost_cell_grid)) {
      cell_part_cnt += cells[c].pList.n;
      cell = &cells[c];
      pl   = &(cell->pList);
      part = pl->part;
      np   = pl->n;
      for(n=0; n<cells[c].pList.n ; n++) {
	if(part[n].r.identity < 0 || part[n].r.identity > max_seen_particle) {
	  fprintf(stderr,"%d: sort_part_in_cells: ERROR: Cell %d Part %d has corrupted id=%d\n",
		  this_node,c,n,cells[c].pList.part[n].r.identity);
	  errexit();
	}
	for(dir=0;dir<3;dir++) {
	  if(periodic[dir] && (part[n].r.p[dir] < 0 || part[n].r.p[dir] > box_l[dir])) {
	    fprintf(stderr,"%d: exchange_part: ERROR: illegal pos[%d]=%f of part %d id=%d in cell %d\n",
		    this_node,dir,part[n].r.p[dir],n,part[n].r.identity,c);
	    errexit();
	  }
	}
	if(local_particles[part[n].r.identity] != &part[n]) {
	    fprintf(stderr,"%d: exchange_part: ERROR: address mismatch for part id %d: %p %p in cell %d\n",
		    this_node,part[n].r.identity,local_particles[part[n].r.identity],
		    &part[n],c);
	    errexit();

	}
      }
    }
    else {
      if(cells[c].pList.n>0) {
	ghost_part_cnt += cells[c].pList.n;
	fprintf(stderr,"%d: sort_part_in_cells: WARNING: ghost_cell %d contains %d particles!\n",
		this_node,c,cells[c].pList.n);
      }
    }
  }
  CELL_TRACE(fprintf(stderr,"%d: sort_part_in_cells: %d particles in cells.\n",
		     this_node,cell_part_cnt));
  /* checks: local particle id */
  for(n=0; n< max_seen_particle+1; n++) {
    if(local_particles[n] != NULL) {
      local_part_cnt ++;
      if(local_particles[n]->r.identity != n) {
	fprintf(stderr,"%d: sort_part_in_cells: ERROR: local_particles part %d has corrupted id %d\n",
		this_node,n,local_particles[n]->r.identity);
	errexit();
      }
    }
  }
  CELL_TRACE(fprintf(stderr,"%d: sort_part_in_cells: %d particles in local_particles.\n",
		     this_node,local_part_cnt));

  /* check cell neighbor consistency */
  for(c=0; c<n_cells; c++) {
    for(n=0; n< cells[c].n_neighbors; n++) {
      if(is_inner_cell(c,ghost_cell_grid)) {
	nc = cells[c].nList[n].cell_ind;
	if( cells[c].nList[n].pList != &(cells[nc].pList) ) {
	  fprintf(stderr,"%d: cell %d: neighbor_cell %d with c_ind %d: Location of pList changed wothout update!\n",this_node,c,n,nc);
	  cell_err_cnt++;
	}
      }
      else {
	fprintf(stderr,"%d: ghost cell %d has more than zero neighbors = %d\n",
		this_node,c,cells[c].n_neighbors);
	cell_err_cnt++;
      }
    }
  }

  /* EXIT on severe errors */
  if(cell_err_cnt>0) {
    fprintf(stderr,"%d: %d ERRORS detected in cell structure!\n",this_node,cell_err_cnt);
    errexit();
  }
  if(local_part_cnt != cell_part_cnt) {
    fprintf(stderr,"%d: sort_part_in_cells: ERROR: %d parts in cells but %d parts in local_particles\n",
	    this_node,local_part_cnt,cell_part_cnt);
    if(ghost_part_cnt==0) errexit();
  }
  if(ghost_part_cnt>0) {
    fprintf(stderr,"%d: sort_part_in_cells: ERROR: Found %d illegal ghost particles!\n",
	    this_node,ghost_part_cnt);
    errexit();
  }
}

