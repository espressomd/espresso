// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file ghosts.c   Ghost particles and particle exchange.
 *
 *  For more information on ghosts,
 *  see \ref ghosts.h "ghosts.h" 
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ghosts.h"
#include "debug.h"
#include "global.h"
#include "cells.h"
#include "communication.h"
#include "grid.h"
#include "utils.h"
#include "particle_data.h"
#include "forces.h"

/** Tag for communication in ghost_comm. */
#define REQ_GHOST_SEND 100

static int n_s_buffer = 0;
static int max_s_buffer = 0;
/** send buffer. Just grows, which should be ok */
static void *s_buffer = NULL;

static int n_r_buffer = 0;
static int max_r_buffer = 0;
/** recv buffer. Just grows, which should be ok */
static void *r_buffer = NULL;

static MPI_Op MPI_FORCES_SUM;

/************************************************************
 * Exported Functions
 ************************************************************/

void prepare_comm(GhostCommunicator *comm, int data_parts, int num)
{
  int i;
  comm->data_parts = data_parts;
  comm->num = num;
  comm->comm = malloc(num*sizeof(GhostCommunication));
  for(i=0; i<num; i++) {
    comm->comm[i].shift[0]=comm->comm[i].shift[1]=comm->comm[i].shift[2]=0.0;
  }
}

void free_comm(GhostCommunicator *comm)
{
  int n;
  GHOST_TRACE(fprintf(stderr,"%d: free_comm: %p has %d ghost communications\n",this_node,comm,comm->num));
  for (n = 0; n < comm->num; n++) free(comm->comm[n].part_lists);
  free(comm->comm);
}

int calc_transmit_size(GhostCommunication *gc, int data_parts)
{
  int p, n_buffer_new;

  if (data_parts == GHOSTTRANS_PARTNUM)
    n_buffer_new = sizeof(int)*gc->n_part_lists;
  else {
    int count = 0;
    for (p = 0; p < gc->n_part_lists; p++)
      count += gc->part_lists[p]->n;

    n_buffer_new = 0;
    if (data_parts & GHOSTTRANS_PROPRTS)
      n_buffer_new += sizeof(ParticleProperties);
    if (data_parts & GHOSTTRANS_POSITION)
      n_buffer_new += sizeof(ParticlePosition);
    if (data_parts & GHOSTTRANS_MOMENTUM)
      n_buffer_new += sizeof(ParticleMomentum);
    if (data_parts & GHOSTTRANS_FORCE)
      n_buffer_new += sizeof(ParticleForce);
    n_buffer_new *= count;
  }
  return n_buffer_new;
}

void prepare_send_buffer(GhostCommunication *gc, int data_parts)
{
  void *insert;
  int pl, p, np;
  Particle *part, *pt;

  /* reallocate send buffer */
  n_s_buffer = calc_transmit_size(gc, data_parts);
  if (n_s_buffer > max_s_buffer) {
    max_s_buffer = n_s_buffer;
    s_buffer = realloc(s_buffer, max_s_buffer);
  }
  /* put in data */
  insert = s_buffer;
  for (pl = 0; pl < gc->n_part_lists; pl++) {
    np   = gc->part_lists[pl]->n;
    if (data_parts == GHOSTTRANS_PARTNUM) {
      *(int *)insert = np;
      insert += sizeof(int);
    }
    else {
      part = gc->part_lists[pl]->part;
      for (p = 0; p < np; p++) {
	pt = &part[p];
	if (data_parts & GHOSTTRANS_PROPRTS) {
	  memcpy(insert, &pt->p, sizeof(ParticleProperties));
	  insert +=  sizeof(ParticleProperties);
	}
	if (data_parts & GHOSTTRANS_POSSHFTD) {
	  /* ok, this is not nice, but perhaps fast */
	  ParticlePosition *pp = insert;
	  int i;
	  memcpy(pp, &pt->r, sizeof(ParticlePosition));
	  for (i = 0; i < 3; i++)
	    pp->p[i] += gc->shift[i];
	  insert +=  sizeof(ParticlePosition);
	}
	else if (data_parts & GHOSTTRANS_POSITION) {
	  memcpy(insert, &pt->r, sizeof(ParticlePosition));
	  insert +=  sizeof(ParticlePosition);
	}
	if (data_parts & GHOSTTRANS_MOMENTUM) {
	  memcpy(insert, &pt->m, sizeof(ParticleMomentum));
	  insert +=  sizeof(ParticleMomentum);
	}
	if (data_parts & GHOSTTRANS_FORCE) {
	  memcpy(insert, &pt->f, sizeof(ParticleForce));
	  insert +=  sizeof(ParticleForce);
	}
      }
    }
  }
#ifdef ADDITIONAL_CHECKS
  if (insert - s_buffer != n_s_buffer) {
    fprintf(stderr, "%d: send buffer size %d differs from what I put in %d\n", this_node, n_s_buffer, insert - s_buffer);
    errexit();
  }
#endif
}

void prepare_recv_buffer(GhostCommunication *gc, int data_parts)
{
  /* reallocate recv buffer */
  n_r_buffer = calc_transmit_size(gc, data_parts);
  if (n_r_buffer > max_r_buffer) {
    max_r_buffer = n_r_buffer;
    r_buffer = realloc(r_buffer, max_r_buffer);
  }
}

void put_recv_buffer(GhostCommunication *gc, int data_parts)
{
  int pl, p, np;
  Particle *part, *pt;
  void *retrieve;

  /* put back data */
  retrieve = r_buffer;
  for (pl = 0; pl < gc->n_part_lists; pl++) {
    np   = gc->part_lists[pl]->n;
    if (data_parts == GHOSTTRANS_PARTNUM) {
      *(int *)retrieve = np;
      retrieve += sizeof(int);
    }
    else {
      part = gc->part_lists[pl]->part;
      for (p = 0; p < np; p++) {
	pt = &part[p];
	if (data_parts & GHOSTTRANS_PROPRTS) {
	  memcpy(retrieve, &pt->p, sizeof(ParticleProperties));
	  retrieve +=  sizeof(ParticleProperties);
	  if (local_particles[pt->p.identity] == NULL)
	    local_particles[pt->p.identity] = pt;
	}
	if (data_parts & GHOSTTRANS_POSITION) {
	  memcpy(retrieve, &pt->r, sizeof(ParticlePosition));
	  retrieve +=  sizeof(ParticlePosition);
	}
	if (data_parts & GHOSTTRANS_MOMENTUM) {
	  memcpy(retrieve, &pt->m, sizeof(ParticleMomentum));
	  retrieve +=  sizeof(ParticleMomentum);
	}
	if (data_parts & GHOSTTRANS_FORCE) {
	  memcpy(retrieve, &pt->f, sizeof(ParticleForce));
	  retrieve +=  sizeof(ParticleForce);
	}
      }
    }
  }
#ifdef ADDITIONAL_CHECKS
  if (retrieve - r_buffer != n_r_buffer) {
    fprintf(stderr, "%d: recv buffer size %d differs from what I put in %d\n", this_node, n_r_buffer, retrieve - r_buffer);
    errexit();
  }
#endif
}

void add_forces_from_recv_buffer(GhostCommunication *gc)
{
  int pl, p, np;
  Particle *part, *pt;
  void *retrieve;

  /* put back data */
  retrieve = r_buffer;
  for (pl = 0; pl < gc->n_part_lists; pl++) {
    np   = gc->part_lists[pl]->n;
    part = gc->part_lists[pl]->part;
    for (p = 0; p < np; p++) {
      pt = &part[p];
      add_force(&pt->f, (ParticleForce *)retrieve);
      retrieve +=  sizeof(ParticleForce);
    }
  }
#ifdef ADDITIONAL_CHECKS
  if (retrieve - r_buffer != n_r_buffer) {
    fprintf(stderr, "%d: recv buffer size %d differs from what I put in %d\n", this_node, n_r_buffer, retrieve - r_buffer);
    errexit();
  }
#endif
}

void reduce_forces_sum(void *add, void *to, int *len, MPI_Datatype *type)
{
  void *cur = add;

#ifdef ADDITIONAL_CHECKS
  if (*type != MPI_BYTE) {
    fprintf(stderr, "%d: transfer data type wrong\n", this_node);
    errexit();
  }
#endif

  for (; cur - add < *len; to++, add++)
    add_force((ParticleForce *)to, (ParticleForce *)add);
}

void ghost_communicator(GhostCommunicator *gc)
{
  int n, n2;
  for (n = 0; n < gc->num; n++) {
    int comm_type = gc->comm[n].type & GHOST_JOBMASK;

    /* prepare send buffer if necessary */
    if ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
	(comm_type == GHOST_BCST && gc->comm[n].node == this_node)) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!(gc->comm[n].type & GHOST_PREFETCH))
	prepare_send_buffer(&gc->comm[n], gc->data_parts);
    }
    else {
      /* we do not send this time, let's look for a prefetch */
      if (gc->comm[n].type & GHOST_PREFETCH) {
	for (n2 = n+1; n2 < gc->num; n2++) {
	  int comm_type2 = gc->comm[n2].type & GHOST_JOBMASK;
	  /* find next action where we send and which has PREFETCH set */
	  if (((comm_type2 == GHOST_SEND) || (comm_type2 == GHOST_RDCE) ||
	       (comm_type2 == GHOST_BCST && gc->comm[n2].node == this_node)) &&
	      gc->comm[n2].type & GHOST_PREFETCH)
	    prepare_send_buffer(&gc->comm[n2], gc->data_parts);
	}
      }
    }

    /* recv buffer for recv and multinode operations to this node */
    if ((comm_type == GHOST_RECV) ||
	(comm_type == GHOST_RDCE && gc->comm[n].node == this_node) ||
	(comm_type == GHOST_BCST))
      prepare_recv_buffer(&gc->comm[n], gc->data_parts);

    /* transfer data */

    /* write back data */
    if (comm_type == GHOST_RECV || comm_type == GHOST_BCST) {
      /* forces have to be added, the rest overwritten */
      if (gc->data_parts == GHOSTTRANS_FORCE)
	add_forces_from_recv_buffer(&gc->comm[n]);
      else
	put_recv_buffer(&gc->comm[n], gc->data_parts);
    }
    else if (comm_type == GHOST_RDCE && gc->comm[n].node == this_node)
      /* exception: for reduce, the addition is integral part of the reduce operation */
      put_recv_buffer(&gc->comm[n], gc->data_parts);
  }
}

void ghost_init()
{
  MPI_Op_create(reduce_forces_sum, 1, &MPI_FORCES_SUM);
}

/** Go through \ref ghost_cells and remove the ghost entries from \ref
    local_particles. Part of \ref dd_exchange_and_sort_particles.*/
void invalidate_ghosts()
{
  Particle *part;
  int c, np, p;
  /* remove ghosts, but keep Real Particles */
  for(c=0; c<ghost_cells.n; c++) {
    part = ghost_cells.cell[c]->part;
    np   = ghost_cells.cell[c]->n;
    for(p=0 ; p<np; p++) {
      /* Particle is stored as ghost in the local_particles array,
	 if the pointer stored there belongs to a ghost celll
	 particle array. */
      if( &(part[p]) == local_particles[part[p].p.identity] ) 
	local_particles[part[p].p.identity] = NULL;
    }
    ghost_cells.cell[c]->n = 0;
  }
}
