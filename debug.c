// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file debug.c
    Implements the malloc replacements as described in \ref debug.h "debug.h". */

#include <mpi.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "communication.h"
#include "debug.h"
#include "cells.h"
#include "grid.h"
#include "integrate.h"

#if defined FORCE_CORE || defined MPI_CORE
int regular_exit = 0;
#else
int regular_exit = 1;
#endif
static int core_done = 0;

#ifdef ONEPART_DEBUG
int check_id =  ONEPART_DEBUG ;
#endif

#ifdef MEM_DEBUG

#undef realloc
#undef malloc
#undef free

void *__realloc(void *old, unsigned int size, char *where, int line)
{
  void *ret;
  if (size == 0) {
    fprintf(stderr, "%d: free %p at %s:%d\n", this_node, old, where, line);
    free(old);
    return NULL;
  }

  ret = realloc(old, size);
  fprintf(stderr, "%d: realloc %p -> %p size %d at %s:%d\n", this_node, old, ret, size, where, line);
  return ret;
}

void *__malloc(unsigned int size, char *where, int line)
{
  void *ret;
  if (size > 0)
    ret = malloc(size);
  else
    ret = NULL;
  fprintf(stderr, "%d: alloc %d -> %p at %s:%d\n", this_node, size, ret, where, line);
  return ret;
}

void __free(void *p, char *where, int line)
{
  fprintf(stderr, "%d: free %p at %s:%d\n", this_node, p, where, line);
  free(p);
}

#endif

int debug_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Debug status { ", (char *) NULL);
#ifdef COMM_DEBUG
  Tcl_AppendResult(interp, "COMM_DEBUG ", (char *) NULL);
#endif
#ifdef INTEG_DEBUG
  Tcl_AppendResult(interp, "INTEG_DEBUG ", (char *) NULL);
#endif
#ifdef CELL_DEBUG
  Tcl_AppendResult(interp, "CELL_DEBUG ", (char *) NULL);
#endif
#ifdef GHOST_DEBUG
  Tcl_AppendResult(interp, "GHOST_DEBUG ", (char *) NULL);
#endif
#ifdef GRID_DEBUG 
  Tcl_AppendResult(interp, "GRID_DEBUG ", (char *) NULL);
#endif
#ifdef VERLET_DEBUG
  Tcl_AppendResult(interp, "VERLET_DEBUG ", (char *) NULL);
#endif
#ifdef PARTICLE_DEBUG
  Tcl_AppendResult(interp, "PARTICLE_DEBUG ", (char *) NULL);
#endif
#ifdef P3M_DEBUG
  Tcl_AppendResult(interp, "P3M_DEBUG ", (char *) NULL);
#endif
#ifdef FFT_DEBUG
  Tcl_AppendResult(interp, "FFT_DEBUG ", (char *) NULL);
#endif
#ifdef RANDOM_DEBUG
  Tcl_AppendResult(interp, "RANDOM_DEBUG ", (char *) NULL);
#endif
#ifdef FORCE_DEBUG
  Tcl_AppendResult(interp, "FORCE_DEBUG ", (char *) NULL);
#endif
#ifdef THERMO_DEBUG
  Tcl_AppendResult(interp, "THERMO_DEBUG ", (char *) NULL);
#endif
#ifdef LJ_DEBUG
  Tcl_AppendResult(interp, "LJ_DEBUG ", (char *) NULL);
#endif
#ifdef ESR_DEBUG
  Tcl_AppendResult(interp, "ESR_DEBUG ", (char *) NULL);
#endif
#ifdef ESK_DEBUG
  Tcl_AppendResult(interp, "ESK_DEBUG ", (char *) NULL);
#endif
#ifdef FENE_DEBUG
  Tcl_AppendResult(interp, "FENE_DEBUG ", (char *) NULL);
#endif
#ifdef GHOST_FORCE_DEBUG
  Tcl_AppendResult(interp, "GHOST_FORCE_DEBUG ", (char *) NULL);
#endif
#ifdef ONEPART_DEBUG
  Tcl_AppendResult(interp, "ONEPART_DEBUG ", (char *) NULL);
#endif
#ifdef STAT_DEBUG
  Tcl_AppendResult(interp, "STAT_DEBUG ", (char *) NULL);
#endif
#ifdef POLY_DEBUG
  Tcl_AppendResult(interp, "POLY_DEBUG ", (char *) NULL);
#endif
#ifdef MPI_CORE
  Tcl_AppendResult(interp, "MPI_CORE ", (char *) NULL);
#endif
#ifdef FORCE_CORE
  Tcl_AppendResult(interp, "FORCE_CORE ", (char *) NULL);
#endif
#ifdef ADDITIONAL_CHECKS
  Tcl_AppendResult(interp, "ADDITIONAL_CHECKS ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "} }", (char *) NULL);
  return (TCL_OK);
}

void core()
{
  if (!core_done && !regular_exit) {
    fprintf(stderr, "%d: forcing core dump on irregular exit (%d / %d) \n", this_node, core_done, regular_exit);
    *(int *)0 = 0;
    /* doesn't work on AMD64 */
    /* kill(getpid(), SIGSEGV); */
    core_done = 1;
  }
}

void check_particle_consistency()
{
  Particle *part;
  Cell *cell;
  int n, np, dir, c, p;
  int cell_part_cnt=0, ghost_part_cnt=0, local_part_cnt=0;
  int cell_err_cnt=0;

  /* checks: part_id, part_pos, local_particles id */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    cell_part_cnt += cell->n;
    part = cell->part;
    np   = cell->n;
    for(n=0; n<cell->n ; n++) {
      if(part[n].p.identity < 0 || part[n].p.identity > max_seen_particle) {
	fprintf(stderr,"%d: check_particle_consistency: ERROR: Cell %d Part %d has corrupted id=%d\n",
		this_node,c,n,cell->part[n].p.identity);
	errexit();
      }
      for(dir=0;dir<3;dir++) {
	if(PERIODIC(dir) && (part[n].r.p[dir] < 0 || part[n].r.p[dir] > box_l[dir])) {
	  fprintf(stderr,"%d: check_particle_consistency: ERROR: illegal pos[%d]=%f of part %d id=%d in cell %d\n",
		  this_node,dir,part[n].r.p[dir],n,part[n].p.identity,c);
	  errexit();
	}
      }
      if(local_particles[part[n].p.identity] != &part[n]) {
	fprintf(stderr,"%d: check_particle_consistency: ERROR: address mismatch for part id %d: local: %p cell: %p in cell %d\n",
		this_node,part[n].p.identity,local_particles[part[n].p.identity],
		&part[n],c);
	errexit();
	
      }
    }
  }

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    if(cell->n>0) {
      ghost_part_cnt += cell->n;
      fprintf(stderr,"%d: check_particle_consistency: WARNING: ghost_cell %d contains %d particles!\n",
	      this_node,c,cell->n);
    }
  }
  CELL_TRACE(fprintf(stderr,"%d: check_particle_consistency: %d particles in cells, %d particles in ghost_cells.\n",
		     this_node,cell_part_cnt, ghost_part_cnt));
  /* checks: local particle id */
  for(n=0; n< max_seen_particle+1; n++) {
    if(local_particles[n] != NULL) {
      local_part_cnt ++;
      if(local_particles[n]->p.identity != n) {
	fprintf(stderr,"%d: check_particle_consistency: ERROR: local_particles part %d has corrupted id %d\n",
		this_node,n,local_particles[n]->p.identity);
	errexit();
      }
    }
  }
  CELL_TRACE(fprintf(stderr,"%d: check_particle_consistency: %d particles in local_particles.\n",
		     this_node,local_part_cnt));

  /* EXIT on severe errors */
  if(cell_err_cnt>0) {
    fprintf(stderr,"%d: check_particle_consistency: %d ERRORS detected in cell structure!\n",this_node,cell_err_cnt);
    errexit();
  }
  if(local_part_cnt != cell_part_cnt) {
    fprintf(stderr,"%d: check_particle_consistency: ERROR: %d parts in cells but %d parts in local_particles\n",
	    this_node,cell_part_cnt,local_part_cnt);

    for (c = 0; c < local_cells.n; c++) {
      for(p = 0; p < local_cells.cell[c]->n; p++)
	fprintf(stderr, "%d: got particle %d in cell %d\n", this_node, local_cells.cell[c]->part[p].p.identity, c);
    }
    
    for(p = 0; p < n_total_particles; p++)
      if (local_particles[p])
	fprintf(stderr, "%d: got particle %d in local_particles\n", this_node, p);

    if(ghost_part_cnt==0) errexit();
  }
  if(ghost_part_cnt>0) {
    fprintf(stderr,"%d: check_particle_consistency: ERROR: Found %d illegal ghost particles!\n",
	    this_node,ghost_part_cnt);
    errexit();
  }
}

void check_particles()
{
  Particle *part;
  int *is_here;
  Cell *cell;
  int n, np, dir, c, p;
  int cell_part_cnt=0, local_part_cnt=0;
  int cell_err_cnt=0;
  double skin2 = (skin != -1) ? skin/2 : 0;

  CELL_TRACE(fprintf(stderr, "%d: entering check_particles\n", this_node));

  /* check the consistency of particle_nodes */
  /* to this aim the array is broadcasted temporarily */
  if (this_node != 0)
    particle_node = malloc((max_seen_particle + 1)*sizeof(int));
  is_here = malloc((max_seen_particle + 1)*sizeof(int));
  memset(is_here, 0, (max_seen_particle + 1)*sizeof(int));

  MPI_Bcast(particle_node, max_seen_particle + 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* checks: part_id, part_pos, local_particles id */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    cell_part_cnt += cell->n;
    part = cell->part;
    np   = cell->n;
    for(n=0; n<cell->n ; n++) {
      if(part[n].p.identity < 0 || part[n].p.identity > max_seen_particle) {
	fprintf(stderr,"%d: check_particles: ERROR: Cell %d Part %d has corrupted id=%d\n",
		this_node,c,n,cell->part[n].p.identity);
	errexit();
      }

      is_here[part[n].p.identity] = 1;

      for(dir=0;dir<3;dir++) {
	if(PERIODIC(dir) && (part[n].r.p[dir] < -skin2 || part[n].r.p[dir] > box_l[dir] + skin2)) {
	  fprintf(stderr,"%d: check_particles: ERROR: illegal pos[%d]=%f of part %d id=%d in cell %d\n",
		  this_node,dir,part[n].r.p[dir],n,part[n].p.identity,c);
	  errexit();
	}
      }
      if(local_particles[part[n].p.identity] != &part[n]) {
	fprintf(stderr,"%d: check_particles: ERROR: address mismatch for part id %d: local: %p cell: %p in cell %d\n",
		this_node,part[n].p.identity,local_particles[part[n].p.identity],
		&part[n],c);
	errexit();
      }
      if (particle_node[part[n].p.identity] != this_node) {
	fprintf(stderr,"%d: check_particles: ERROR: node for particle %d wrong\n",
		this_node,part[n].p.identity);
	errexit();
      }
    }
  }
  CELL_TRACE(fprintf(stderr,"%d: check_particles: %d particles in local cells.\n",
		     this_node,cell_part_cnt));

  /* checks: local particle id */
  for(n = 0; n <= max_seen_particle; n++) {
    if(local_particles[n] != NULL) {
      local_part_cnt ++;
      if(local_particles[n]->p.identity != n) {
	fprintf(stderr,"%d: check_particles: ERROR: local_particles part %d has corrupted id %d\n",
		this_node,n,local_particles[n]->p.identity);
	errexit();
      }
    }
  }
  CELL_TRACE(fprintf(stderr,"%d: check_particles: %d particles in local_particles.\n",
		     this_node,local_part_cnt));

  /* EXIT on severe errors */
  if(cell_err_cnt>0) {
    fprintf(stderr,"%d: check_particles: %d ERRORS detected in cell structure!\n",this_node,cell_err_cnt);
    errexit();
  }

  /* check whether the particles on my node are actually here */
  for (p = 0; p <= max_seen_particle; p++) {
    if (particle_node[p] == this_node) {
      if (!is_here[p]) {
 	fprintf(stderr,"%d: check_particles: ERROR: particle %d on this node, but not in local cell\n", this_node, p);
      }
    }
  }

  free(is_here);

  if (this_node != 0) {
    free(particle_node);
    particle_node = NULL;
  }
  else {
    /* check whether the total count of particles is ok */
    c = 0;
    for (p = 0; p <= max_seen_particle; p++)
      if (particle_node[p] != -1) c++;
    if (c != n_total_particles) {
      fprintf(stderr,"%d: check_particles: #particles in particle_node inconsistent\n", this_node);
      errexit();
    }
    CELL_TRACE(fprintf(stderr,"%d: check_particles: %d particles in particle_node.\n",
		       this_node,c));
  }
  CELL_TRACE(fprintf(stderr, "%d: leaving check_particles\n", this_node));
}

void print_particle_positions()
{
  int c,np,p;
  Cell *cell;
  Particle *part;

  for(c=0;c<n_cells;c++) {
    cell = &cells[c];
    part = cell->part;
    np = cell->n;
    if(np>0) {
      fprintf(stderr,"%d: cell %d contains %d particles:\n",this_node,c,np);
      for(p=0;p<np;p++) {
	fprintf(stderr,"%d: c%d p%d (id%d) pos %f %f %f\n",this_node,c,p,part[p].p.identity,part[p].r.p[0],part[p].r.p[1],part[p].r.p[2]);
      }
    }
  }    
}

void print_particle_forces()
{
  int c,np,p;
  Cell *cell;
  Particle *part;

  for(c=0;c<n_cells;c++) {
    cell = &cells[c];
    part = cell->part;
    np = cell->n;
    if(np>0) {
      fprintf(stderr,"%d: cell %d contains %d particles:\n",this_node,c,np);
      for(p=0;p<np;p++) {
	fprintf(stderr,"%d: c%d p%d (id%d) force %f %f %f\n",this_node,c,p,part[p].p.identity,part[p].f.f[0],part[p].f.f[1],part[p].f.f[2]);
      }
    }
  }    
}
