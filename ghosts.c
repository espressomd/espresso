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
static char *s_buffer = NULL;

static int n_r_buffer = 0;
static int max_r_buffer = 0;
/** recv buffer. Just grows, which should be ok */
static char *r_buffer = NULL;

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
  char *insert;
  int pl, p, np;
  Particle *part, *pt;

  GHOST_TRACE(fprintf(stderr, "%d: prepare sending to/bcast from %d\n", this_node, gc->node));

  /* reallocate send buffer */
  n_s_buffer = calc_transmit_size(gc, data_parts);
  if (n_s_buffer > max_s_buffer) {
    max_s_buffer = n_s_buffer;
    s_buffer = realloc(s_buffer, max_s_buffer);
  }
  GHOST_TRACE(fprintf(stderr, "%d: will send %d\n", this_node, n_s_buffer));

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
	  ParticlePosition *pp = (ParticlePosition *)insert;
	  int i;
	  memcpy(pp, &pt->r, sizeof(ParticlePosition));
	  //fprintf(stderr,"%d prep_send_buf: ghost %d shift %f,%f,%f\n",this_node,pt->p.identity,gc->shift[0],gc->shift[1],gc->shift[2]);
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
  GHOST_TRACE(fprintf(stderr, "%d: prepare receiving from %d\n", this_node, gc->node));
  /* reallocate recv buffer */
  n_r_buffer = calc_transmit_size(gc, data_parts);
  if (n_r_buffer > max_r_buffer) {
    max_r_buffer = n_r_buffer;
    r_buffer = realloc(r_buffer, max_r_buffer);
  }
  GHOST_TRACE(fprintf(stderr, "%d: will get %d\n", this_node, n_r_buffer));
}

void put_recv_buffer(GhostCommunication *gc, int data_parts)
{
  int pl, p, np;
  Particle *part, *pt;
  char *retrieve;

  /* put back data */
  retrieve = r_buffer;
  for (pl = 0; pl < gc->n_part_lists; pl++) {
    if (data_parts == GHOSTTRANS_PARTNUM) {
      GHOST_TRACE(fprintf(stderr, "%d: reallocating cell %p to size %d\n", this_node, gc->part_lists[pl], *(int *)retrieve));
      realloc_particlelist(gc->part_lists[pl], gc->part_lists[pl]->n = *(int *)retrieve);
      retrieve += sizeof(int);
    }
    else {
      np   = gc->part_lists[pl]->n;
      part = gc->part_lists[pl]->part;
      for (p = 0; p < np; p++) {
	pt = &part[p];
	if (data_parts & GHOSTTRANS_PROPRTS) {
	  memcpy(&pt->p, retrieve, sizeof(ParticleProperties));
	  retrieve +=  sizeof(ParticleProperties);
	  /* GHOST_TRACE(fprintf(stderr, "%d: received ghost %d", this_node, pt->p.identity)); */
	  if (local_particles[pt->p.identity] == NULL) {
	    /* GHOST_TRACE(fprintf(stderr, ", using.\n")); */
	    local_particles[pt->p.identity] = pt;
	  }
	  /*
	    else {
	    GHOST_TRACE(fprintf(stderr, ", already here.\n"));
	    }
	  */
	}
	if (data_parts & GHOSTTRANS_POSITION) {
	  memcpy(&pt->r, retrieve, sizeof(ParticlePosition));
	  retrieve +=  sizeof(ParticlePosition);
	}
	if (data_parts & GHOSTTRANS_MOMENTUM) {
	  memcpy(&pt->m, retrieve, sizeof(ParticleMomentum));
	  retrieve +=  sizeof(ParticleMomentum);
	}
	if (data_parts & GHOSTTRANS_FORCE) {
	  memcpy(&pt->f, retrieve, sizeof(ParticleForce));
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
  char *retrieve;

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

void cell_cell_transfer(GhostCommunication *gc, int data_parts)
{
  int pl, p, np, offset;
  Particle *part1, *part2, *pt1, *pt2;

  GHOST_TRACE(fprintf(stderr, "%d: local_transfer: type %d data_parts %d\n", this_node,gc->type,data_parts));

  /* transfer data */
  offset = gc->n_part_lists/2;
  for (pl = 0; pl < offset; pl++) {
    np   = gc->part_lists[pl]->n;
    if (data_parts == GHOSTTRANS_PARTNUM) {
      realloc_particlelist(gc->part_lists[pl + offset],
			   gc->part_lists[pl + offset]->n = gc->part_lists[pl]->n); 
    }
    else {
      part1 = gc->part_lists[pl]->part;
      part2 = gc->part_lists[pl + offset]->part;
      for (p = 0; p < np; p++) {
	pt1 = &part1[p];
	pt2 = &part2[p];
	if (data_parts & GHOSTTRANS_PROPRTS)
	  memcpy(&pt2->p, &pt1->p, sizeof(ParticleProperties));
	if (data_parts & GHOSTTRANS_POSSHFTD) {
	  /* ok, this is not nice, but perhaps fast */
	  int i;
	  memcpy(&pt2->r, &pt1->r, sizeof(ParticlePosition));
	  for (i = 0; i < 3; i++)
	    pt2->r.p[i] += gc->shift[i];
	}
	else if (data_parts & GHOSTTRANS_POSITION)
	  memcpy(&pt2->r, &pt1->r, sizeof(ParticlePosition));
	if (data_parts & GHOSTTRANS_MOMENTUM)
	  memcpy(&pt2->m, &pt1->m, sizeof(ParticleMomentum));
	if (data_parts & GHOSTTRANS_FORCE)
	  add_force(&pt2->f, &pt1->f);
      }
    }
  }
}

void reduce_forces_sum(void *add, void *to, int *len, MPI_Datatype *type)
{
  ParticleForce *cadd = add, *cto = to;
  int i, clen = *len/sizeof(ParticleForce);
 
#ifdef ADDITIONAL_CHECKS
  if (*type != MPI_BYTE || (*len % sizeof(ParticleForce)) != 0) {
    fprintf(stderr, "%d: transfer data type wrong\n", this_node);
    errexit();
  }
#endif

  for (i = 0; i < clen; i++)
    add_force(&cto[i], &cadd[i]);
}

static int is_send_op(int comm_type, int node)
{
  return ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
	  (comm_type == GHOST_BCST && node == this_node));
}

static int is_recv_op(int comm_type, int node)
{
  return ((comm_type == GHOST_RECV) ||
	  (comm_type == GHOST_BCST && node != this_node) ||
	  (comm_type == GHOST_RDCE && node == this_node));
}

void ghost_communicator(GhostCommunicator *gc)
{
  MPI_Status status;
  int n, n2;
  int data_parts = gc->data_parts;

  GHOST_TRACE(fprintf(stderr, "%d: ghost_comm %p, data_parts %d\n", this_node, gc, data_parts));

  for (n = 0; n < gc->num; n++) {
    GhostCommunication *gcn = &gc->comm[n];
    int comm_type = gcn->type & GHOST_JOBMASK;
    int prefetch  = gcn->type & GHOST_PREFETCH;
    int node      = gcn->node;
    
    GHOST_TRACE(fprintf(stderr, "%d: round %d, job %x\n", this_node, n, gc->comm[n].type));
    GHOST_TRACE(fprintf(stderr, "%d: shift %f %f %f\n",this_node, gc->comm[n].shift[0], gc->comm[n].shift[1], gc->comm[n].shift[2]));
    if (comm_type == GHOST_LOCL)
      cell_cell_transfer(gcn, data_parts);
    else {
      /* prepare send buffer if necessary */
      if (is_send_op(comm_type, node)) {
	/* ok, we send this step, prepare send buffer if not yet done */
	if (!prefetch)
	  prepare_send_buffer(gcn, data_parts);
      }
      else {
	/* we do not send this time, let's look for a prefetch */
	if (prefetch) {
	  for (n2 = n+1; n2 < gc->num; n2++) {
	    GhostCommunication *gcn2 = &gc->comm[n2];
	    int comm_type2 = gcn2->type & GHOST_JOBMASK;
	    int prefetch2  = gcn2->type & GHOST_PREFETCH;
	    int node2      = gcn2->node;
	    /* find next action where we send and which has PREFETCH set */
	    if (is_send_op(comm_type2, node2) && prefetch2) {
	      GHOST_TRACE(fprintf(stderr, "%d: prefetch comm with/from %d\n", this_node, node2));
	      prepare_send_buffer(gcn2, data_parts);
	    }
	  }
	}
      }

      /* recv buffer for recv and multinode operations to this node */
      if (is_recv_op(comm_type, node))
	prepare_recv_buffer(gcn, data_parts);

      /* transfer data */
      switch (comm_type) {
      case GHOST_RECV:
	GHOST_TRACE(fprintf(stderr, "%d: receive from %d\n", this_node, node));
	MPI_Recv(r_buffer, n_r_buffer, MPI_BYTE, node, REQ_GHOST_SEND, MPI_COMM_WORLD, &status);
	break;
      case GHOST_SEND:
	GHOST_TRACE(fprintf(stderr, "%d: send to %d\n", this_node, node));
	MPI_Send(s_buffer, n_s_buffer, MPI_BYTE, node, REQ_GHOST_SEND, MPI_COMM_WORLD);
	break;
      case GHOST_BCST:
	GHOST_TRACE(fprintf(stderr, "%d: bcast from %d\n", this_node, node));
	if (node == this_node)
	  MPI_Bcast(s_buffer, n_s_buffer, MPI_BYTE, node, MPI_COMM_WORLD);
	else
	  MPI_Bcast(r_buffer, n_r_buffer, MPI_BYTE, node, MPI_COMM_WORLD);
	break;
      case GHOST_RDCE:
	GHOST_TRACE(fprintf(stderr, "%d: reduce to %d\n", this_node, node));
	if (node == this_node)
	  MPI_Reduce(s_buffer, r_buffer, n_s_buffer, MPI_BYTE, MPI_FORCES_SUM, node, MPI_COMM_WORLD);
	else
	  MPI_Reduce(s_buffer, NULL, n_s_buffer, MPI_BYTE, MPI_FORCES_SUM, node, MPI_COMM_WORLD);
	break;
      }
      GHOST_TRACE(fprintf(stderr, "%d: done\n", this_node));

      /* write back data */
      if (comm_type == GHOST_RECV ||
	  (comm_type == GHOST_BCST && node != this_node)) {
	/* forces have to be added, the rest overwritten */
	if (data_parts == GHOSTTRANS_FORCE)
	  add_forces_from_recv_buffer(gcn);
	else
	  put_recv_buffer(gcn, data_parts);
      }
      else if (comm_type == GHOST_RDCE && node == this_node)
	/* exception: for reduce, the addition is integral part of the reduce operation */
	put_recv_buffer(gcn, data_parts);
    }
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
