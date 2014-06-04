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
/** \file ghosts.cpp   Ghost particles and particle exchange.
 *
 *  For more information on ghosts,
 *  see \ref ghosts.hpp "ghosts.h" 
*/
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include "utils.hpp"
#include "ghosts.hpp"
#include "global.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "forces.hpp"

/** Tag for communication in ghost_comm. */
#define REQ_GHOST_SEND 100

static int n_s_buffer = 0;
static int max_s_buffer = 0;
/** send buffer. Just grows, which should be ok */
static char *s_buffer = NULL;

std::vector<int> s_bondbuffer;

static int n_r_buffer = 0;
static int max_r_buffer = 0;
/** recv buffer. Just grows, which should be ok */
static char *r_buffer = NULL;

std::vector<int> r_bondbuffer;

static MPI_Op MPI_FORCES_SUM;

/** wether the ghosts should also have velocity information, e. g. for DPD or RATTLE.
    You need this whenever you need the relative velocity of two particles.
    NO CHANGES OF THIS VALUE OUTSIDE OF \ref on_ghost_flags_change !!!!
*/
int ghosts_have_v = 0;

void prepare_comm(GhostCommunicator *comm, int data_parts, int num)
{
  int i;
  comm->data_parts = data_parts;

  /* if ghosts should have uptodate velocities, they have to be updated like positions
     (except for shifting...) */
  if (ghosts_have_v && (data_parts & GHOSTTRANS_POSITION))
    comm->data_parts |= GHOSTTRANS_MOMENTUM;

  GHOST_TRACE(fprintf(stderr, "%d: prepare_comm, data_parts = %d\n", this_node, comm->data_parts));

  comm->num = num;
  comm->comm = (GhostCommunication*)malloc(num*sizeof(GhostCommunication));
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
    if (data_parts & GHOSTTRANS_PROPRTS) {
      n_buffer_new += sizeof(ParticleProperties);
      // sending size of bond/exclusion lists
#ifdef GHOSTS_HAVE_BONDS
      n_buffer_new += sizeof(int);
#ifdef EXCLUSIONS
      n_buffer_new += sizeof(int);
#endif
#endif
    }
    if (data_parts & GHOSTTRANS_POSITION)
      n_buffer_new += sizeof(ParticlePosition);
    if (data_parts & GHOSTTRANS_MOMENTUM)
      n_buffer_new += sizeof(ParticleMomentum);
    if (data_parts & GHOSTTRANS_FORCE)
      n_buffer_new += sizeof(ParticleForce);
#ifdef LB
    if (data_parts & GHOSTTRANS_COUPLING)
      n_buffer_new += sizeof(ParticleLatticeCoupling);
#endif
    n_buffer_new *= count;
  }
  // also sending length of bond buffer
  if (data_parts & GHOSTTRANS_PROPRTS)
    n_buffer_new += sizeof(int);
  return n_buffer_new;
}

void prepare_send_buffer(GhostCommunication *gc, int data_parts)
{
  GHOST_TRACE(fprintf(stderr, "%d: prepare sending to/bcast from %d\n", this_node, gc->node));

  /* reallocate send buffer */
  n_s_buffer = calc_transmit_size(gc, data_parts);
  if (n_s_buffer > max_s_buffer) {
    max_s_buffer = n_s_buffer;
    s_buffer = (char*)realloc(s_buffer, max_s_buffer);
  }
  GHOST_TRACE(fprintf(stderr, "%d: will send %d\n", this_node, n_s_buffer));

  s_bondbuffer.resize(0);

  /* put in data */
  char *insert = s_buffer;
  for (int pl = 0; pl < gc->n_part_lists; pl++) {
    int np   = gc->part_lists[pl]->n;
    if (data_parts == GHOSTTRANS_PARTNUM) {
      *(int *)insert = np;
      insert += sizeof(int);
      GHOST_TRACE(fprintf(stderr, "%d: %d particles assigned\n",
			  this_node, np));
    }
    else {
      Particle *part = gc->part_lists[pl]->part;
      for (int p = 0; p < np; p++) {
	Particle *pt = &part[p];
	if (data_parts & GHOSTTRANS_PROPRTS) {
	  memcpy(insert, &pt->p, sizeof(ParticleProperties));
	  insert +=  sizeof(ParticleProperties);
#ifdef GHOSTS_HAVE_BONDS
          *(int *)insert = pt->bl.n;
	  insert += sizeof(int);
          if (pt->bl.n) {
            s_bondbuffer.insert(s_bondbuffer.end(), pt->bl.e, pt->bl.e + pt->bl.n);
          }
#ifdef EXCLUSIONS
          *(int *)insert = pt->el.n;
	  insert += sizeof(int);
          if (pt->el.n) {
            s_bondbuffer.insert(s_bondbuffer.end(), pt->el.e, pt->el.e + pt->el.n);
          }
#endif
#endif
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
#ifdef LB
	if (data_parts & GHOSTTRANS_COUPLING) {
	  memcpy(insert, &pt->lc, sizeof(ParticleLatticeCoupling));
	  insert +=  sizeof(ParticleLatticeCoupling);
	}
#endif
      }
    }
  }
  if (data_parts & GHOSTTRANS_PROPRTS) {
    GHOST_TRACE(fprintf(stderr, "%d: bond buffer size is %ld\n",
                        this_node, s_bondbuffer.size()));
    *(int *)insert = int(s_bondbuffer.size());
    insert += sizeof(int);
  }

  if (insert - s_buffer != n_s_buffer) {
    fprintf(stderr, "%d: INTERNAL ERROR: send buffer size %d "
            "differs from what I put in (%ld)\n", 
            this_node, n_s_buffer, insert - s_buffer);
    errexit();
  }
}

static void prepare_ghost_cell(Cell *cell, int size)
{
#ifdef GHOSTS_HAVE_BONDS
  // free all allocated information, will be resent
  {
    int np   = cell->n;
    Particle *part = cell->part;
    for (int p = 0; p < np; p++) {
      free_particle(part + p);
    }
  }          
#endif
  realloc_particlelist(cell, cell->n = size);
  // invalidate pointers etc
  {
    int np   = cell->n;
    Particle *part = cell->part;
    for (int p = 0; p < np; p++) {
      Particle *pt = &part[p];
      // no bonds or exclusions
      pt->bl.e = 0;
      pt->bl.n = 0;
      pt->bl.max = 0;
#ifdef EXCLUSIONS
      pt->el.e = 0;
      pt->el.n = 0;
      pt->el.max = 0;
#endif
#ifdef GHOST_FLAG
      //init ghost variable
      pt->l.ghost=1;
#endif
    }
  }
}

void prepare_recv_buffer(GhostCommunication *gc, int data_parts)
{
  GHOST_TRACE(fprintf(stderr, "%d: prepare receiving from %d\n", this_node, gc->node));
  /* reallocate recv buffer */
  n_r_buffer = calc_transmit_size(gc, data_parts);
  if (n_r_buffer > max_r_buffer) {
    max_r_buffer = n_r_buffer;
    r_buffer = (char*)realloc(r_buffer, max_r_buffer);
  }
  GHOST_TRACE(fprintf(stderr, "%d: will get %d\n", this_node, n_r_buffer));
}

void put_recv_buffer(GhostCommunication *gc, int data_parts)
{
  /* put back data */
  char *retrieve = r_buffer;

  std::vector<int>::const_iterator bond_retrieve = r_bondbuffer.begin();

  for (int pl = 0; pl < gc->n_part_lists; pl++) {
    ParticleList *cur_list = gc->part_lists[pl];
    if (data_parts == GHOSTTRANS_PARTNUM) {
      GHOST_TRACE(fprintf(stderr, "%d: reallocating cell %p to size %d, assigned to node %d\n",
			  this_node, cur_list, *(int *)retrieve, gc->node));
      prepare_ghost_cell(cur_list, *(int *)retrieve);
      retrieve += sizeof(int);
    }
    else {
      int np   = cur_list->n;
      Particle *part = cur_list->part;
      for (int p = 0; p < np; p++) {
	Particle *pt = &part[p];
	if (data_parts & GHOSTTRANS_PROPRTS) {
	  memcpy(&pt->p, retrieve, sizeof(ParticleProperties));
	  retrieve +=  sizeof(ParticleProperties);
#ifdef GHOSTS_HAVE_BONDS
          int n_bonds;
	  memcpy(&n_bonds, retrieve, sizeof(int));
	  retrieve +=  sizeof(int);
          if (n_bonds) {
            realloc_intlist(&pt->bl, pt->bl.n = n_bonds);
            std::copy(bond_retrieve, bond_retrieve + n_bonds, pt->bl.e);
            bond_retrieve += n_bonds;
          }
#ifdef EXCLUSIONS
	  memcpy(&n_bonds, retrieve, sizeof(int));
	  retrieve +=  sizeof(int);
          if (n_bonds) {
            realloc_intlist(&pt->el, pt->el.n = n_bonds);
            std::copy(bond_retrieve, bond_retrieve + n_bonds, pt->el.e);
            bond_retrieve += n_bonds;
          }
#endif
#endif
	  if (local_particles[pt->p.identity] == NULL) {
	    local_particles[pt->p.identity] = pt;
	  }
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
#ifdef LB
	if (data_parts & GHOSTTRANS_COUPLING) {
	  memcpy(&pt->lc, retrieve, sizeof(ParticleLatticeCoupling));
	  retrieve +=  sizeof(ParticleLatticeCoupling);
	}
#endif
      }
    }
  }

  if (data_parts & GHOSTTRANS_PROPRTS) {
    // skip the final information on bonds to be sent in a second round
    retrieve += sizeof(int);
  }

  if (retrieve - r_buffer != n_r_buffer) {
    fprintf(stderr, "%d: recv buffer size %d differs "
            "from what I read out (%ld)\n", 
            this_node, n_r_buffer, retrieve - r_buffer);
    errexit();
  }
  if (bond_retrieve != r_bondbuffer.end()) {
    fprintf(stderr, "%d: recv bond buffer was not used up, %ld elements remain\n", 
            this_node, r_bondbuffer.end() - bond_retrieve );
    errexit();
  }
  r_bondbuffer.resize(0);
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
  if (retrieve - r_buffer != n_r_buffer) {
    fprintf(stderr, "%d: recv buffer size %d differs "
            "from what I put in %ld\n", 
            this_node, n_r_buffer, retrieve - r_buffer);
    errexit();
  }
}

void cell_cell_transfer(GhostCommunication *gc, int data_parts)
{
  int pl, p, offset;
  Particle *part1, *part2, *pt1, *pt2;

  GHOST_TRACE(fprintf(stderr, "%d: local_transfer: type %d data_parts %d\n", this_node,gc->type,data_parts));

  /* transfer data */
  offset = gc->n_part_lists/2;
  for (pl = 0; pl < offset; pl++) {
    Cell *src_list = gc->part_lists[pl];
    Cell *dst_list = gc->part_lists[pl + offset];
    if (data_parts == GHOSTTRANS_PARTNUM) {
      prepare_ghost_cell(dst_list, src_list->n);
    }
    else {
      int np = src_list->n;
      part1 = src_list->part;
      part2 = dst_list->part;
      for (p = 0; p < np; p++) {
	pt1 = &part1[p];
	pt2 = &part2[p];
	if (data_parts & GHOSTTRANS_PROPRTS) {
	  memcpy(&pt2->p, &pt1->p, sizeof(ParticleProperties));
#ifdef GHOSTS_HAVE_BONDS
          realloc_intlist(&(pt2->bl), pt2->bl.n = pt1->bl.n);
	  memcpy(&pt2->bl.e, &pt1->bl.e, pt1->bl.n*sizeof(int));
#ifdef EXCLUSIONS
          realloc_intlist(&(pt2->el), pt2->el.n = pt1->el.n);
	  memcpy(&pt2->el.e, &pt1->el.e, pt1->el.n*sizeof(int));
#endif
#endif
        }
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
#ifdef LB
	if (data_parts & GHOSTTRANS_COUPLING)
	  memcpy(&pt2->lc, &pt1->lc, sizeof(ParticleLatticeCoupling));
#endif
      }
    }
  }
}

void reduce_forces_sum(void *add, void *to, int *len, MPI_Datatype *type)
{
  ParticleForce 
    *cadd = (ParticleForce*)add, 
    *cto = (ParticleForce*)to;
  int i, clen = *len/sizeof(ParticleForce);
 
  if (*type != MPI_BYTE || (*len % sizeof(ParticleForce)) != 0) {
    fprintf(stderr, "%d: transfer data type wrong\n", this_node);
    errexit();
  }

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
    int poststore = gcn->type & GHOST_PSTSTORE;
    int node      = gcn->node;
    
    GHOST_TRACE(fprintf(stderr, "%d: ghost_comm round %d, job %x\n", this_node, n, gc->comm[n].type));
    GHOST_TRACE(fprintf(stderr, "%d: ghost_comm shift %f %f %f\n",this_node, gc->comm[n].shift[0], gc->comm[n].shift[1], gc->comm[n].shift[2]));
    if (comm_type == GHOST_LOCL)
      cell_cell_transfer(gcn, data_parts);
    else {
      /* prepare send buffer if necessary */
      if (is_send_op(comm_type, node)) {
	/* ok, we send this step, prepare send buffer if not yet done */
	if (!prefetch)
	  prepare_send_buffer(gcn, data_parts);
	else {
	  GHOST_TRACE(fprintf(stderr, "%d: ghost_comm using prefetched data for operation %d, sending to %d\n", this_node, n, node));
#ifdef ADDITIONAL_CHECKS
	  if (n_s_buffer != calc_transmit_size(gcn, data_parts)) {
	    fprintf(stderr, "%d: ghost_comm transmission size and current size of cells to transmit do not match\n", this_node);
	    errexit();
	  }
#endif
	}
      }
      else {
	/* we do not send this time, let's look for a prefetch */
	if (prefetch) {
	  /* find next action where we send and which has PREFETCH set */
	  for (n2 = n+1; n2 < gc->num; n2++) {
	    GhostCommunication *gcn2 = &gc->comm[n2];
	    int comm_type2 = gcn2->type & GHOST_JOBMASK;
	    int prefetch2  = gcn2->type & GHOST_PREFETCH;
	    int node2      = gcn2->node;
	    if (is_send_op(comm_type2, node2) && prefetch2) {
	      GHOST_TRACE(fprintf(stderr, "%d: ghost_comm prefetch operation %d, is send/bcast to/from %d\n", this_node, n2, node2));
	      prepare_send_buffer(gcn2, data_parts);
	      break;
	    }
	  }
	}
      }

      /* recv buffer for recv and multinode operations to this node */
      if (is_recv_op(comm_type, node))
	prepare_recv_buffer(gcn, data_parts);

      /* transfer data */
      switch (comm_type) {
      case GHOST_RECV: {
	GHOST_TRACE(fprintf(stderr, "%d: ghost_comm receive from %d (%d bytes)\n", this_node, node, n_r_buffer));
	MPI_Recv(r_buffer, n_r_buffer, MPI_BYTE, node, REQ_GHOST_SEND, comm_cart, &status);
        if (data_parts & GHOSTTRANS_PROPRTS) {
          int n_bonds = *(int *)(r_buffer + n_r_buffer - sizeof(int));
          GHOST_TRACE(fprintf(stderr, "%d: ghost_comm receive from %d (%d bonds)\n", this_node, node, n_bonds));
          if (n_bonds) {
            r_bondbuffer.resize(n_bonds);
            MPI_Recv(&r_bondbuffer[0], n_bonds, MPI_INT, node, REQ_GHOST_SEND, comm_cart, &status);
          }
        }
	break;
      }
      case GHOST_SEND: {
	GHOST_TRACE(fprintf(stderr, "%d: ghost_comm send to %d (%d bytes)\n", this_node, node, n_s_buffer));
	MPI_Send(s_buffer, n_s_buffer, MPI_BYTE, node, REQ_GHOST_SEND, comm_cart);
        int n_bonds = s_bondbuffer.size();
        if (!(data_parts & GHOSTTRANS_PROPRTS) && n_bonds > 0) {
          fprintf(stderr, "%d: INTERNAL ERROR: not sending properties, but bond buffer not empty\n", this_node);
          errexit();
        }
        GHOST_TRACE(fprintf(stderr, "%d: ghost_comm send to %d (%d ints)\n", this_node, node, n_bonds));
        if (n_bonds) {
          MPI_Send(&s_bondbuffer[0], n_bonds, MPI_INT, node, REQ_GHOST_SEND, comm_cart);
        }
        break;
      }
      case GHOST_BCST:
	GHOST_TRACE(fprintf(stderr, "%d: ghost_comm bcast from %d (%d bytes)\n", this_node, node,
			    (node == this_node) ? n_s_buffer : n_r_buffer));
	if (node == this_node) {
	  MPI_Bcast(s_buffer, n_s_buffer, MPI_BYTE, node, comm_cart);
          int n_bonds = s_bondbuffer.size();
          if (!(data_parts & GHOSTTRANS_PROPRTS) && n_bonds > 0) {
            fprintf(stderr, "%d: INTERNAL ERROR: not sending properties, but bond buffer not empty\n", this_node);
            errexit();
          }
          if (n_bonds) {
            MPI_Bcast(&s_bondbuffer[0], n_bonds, MPI_INT, node, comm_cart);
          }
        }
	else {
	  MPI_Bcast(r_buffer, n_r_buffer, MPI_BYTE, node, comm_cart);
          if (data_parts & GHOSTTRANS_PROPRTS) {
            int n_bonds = *(int *)(r_buffer + n_r_buffer - sizeof(int));
            if (n_bonds) {
              r_bondbuffer.resize(n_bonds);
              MPI_Bcast(&r_bondbuffer[0], n_bonds, MPI_INT, node, comm_cart);
            }
          }
        }
	break;
      case GHOST_RDCE:
	GHOST_TRACE(fprintf(stderr, "%d: ghost_comm reduce to %d (%d bytes)\n", this_node, node, n_s_buffer));
	if (node == this_node)
	  MPI_Reduce(s_buffer, r_buffer, n_s_buffer, MPI_BYTE, MPI_FORCES_SUM, node, comm_cart);
	else
	  MPI_Reduce(s_buffer, NULL, n_s_buffer, MPI_BYTE, MPI_FORCES_SUM, node, comm_cart);
	break;
      }
      GHOST_TRACE(MPI_Barrier(comm_cart));
      GHOST_TRACE(fprintf(stderr, "%d: ghost_comm done\n", this_node));

      /* recv op; write back data directly, if no PSTSTORE delay is requested. */
      if (is_recv_op(comm_type, node)) {
	if (!poststore) {
	  /* forces have to be added, the rest overwritten. Exception is RDCE, where the addition
	     is integrated into the communication. */
	  if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
	    add_forces_from_recv_buffer(gcn);
	  else
	    put_recv_buffer(gcn, data_parts);
	}
	else {
	  GHOST_TRACE(fprintf(stderr, "%d: ghost_comm delaying operation %d, recv from %d\n", this_node, n, node));
	}
      }
      else {
	/* send op; write back delayed data from last recv, when this was a prefetch send. */
	if (poststore) {
	  /* find previous action where we recv and which has PSTSTORE set */
	  for (n2 = n-1; n2 >= 0; n2--) {
	    GhostCommunication *gcn2 = &gc->comm[n2];
	    int comm_type2 = gcn2->type & GHOST_JOBMASK;
	    int poststore2 = gcn2->type & GHOST_PSTSTORE;
	    int node2      = gcn2->node;
	    if (is_recv_op(comm_type2, node2) && poststore2) {
	      GHOST_TRACE(fprintf(stderr, "%d: ghost_comm storing delayed recv, operation %d, from %d\n", this_node, n2, node2));
#ifdef ADDITIONAL_CHECKS
	      if (n_r_buffer != calc_transmit_size(gcn2, data_parts)) {
		fprintf(stderr, "%d: ghost_comm transmission size and current size of cells to transmit do not match\n", this_node);
		errexit();
	      }
#endif
	      /* as above */
	      if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
		add_forces_from_recv_buffer(gcn2);
	      else
		put_recv_buffer(gcn2, data_parts);
	      break;
	    }
	  }
	}
      }
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
