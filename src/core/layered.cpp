/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file layered.cpp
    Implementation of \ref layered.hpp "layered.h".
 */
#include <mpi.h>
#include <cstring>
#include "utils.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "layered.hpp"
#include "global.hpp"
#include "communication.hpp"
#include "ghosts.hpp"
#include "forces.hpp"
#include "pressure.hpp"
#include "energy.hpp"
#include "constraint.hpp"
#include "domain_decomposition.hpp"

/* Organization: Layers only in one direction.
   ghost_bottom
   c1
   c2
   c3
   .
   .
   .
   cn
   ghost_top

   First, all nodes send downwards, then upwards. Within these actions,
   first the odd nodes send. For even n_nodes, this algorithm is straight
   forward: first the odd nodes send, the even receive, then vice versa.
   For odd n_nodes, we have
   1) 1->2 3->4 5->1
   2) 2->3 4->5
   So in the first round node 5 has to wait for node 1 to
   complete the send and get ready to receive. In other words,
   what physically happens is:
   1) 1->2 3->4 5->*
   2) *->1 2->3 4->*
   3) *->5
   This means that one pending receive has to be done in addition
   provided that all other transmissions can happen in parallel.

*/

/** wether we are the lowest node */
#define LAYERED_BOTTOM 1
/** wether we are the highest node */
#define LAYERED_TOP    2
/** same as PERIODIC(2) */
#define LAYERED_PERIODIC 4
#define LAYERED_BTM_MASK (LAYERED_BOTTOM|LAYERED_PERIODIC)
#define LAYERED_TOP_MASK (LAYERED_TOP|LAYERED_PERIODIC)
/** node has a neighbor above (modulo n_nodes) */
#define LAYERED_TOP_NEIGHBOR ((layered_flags & LAYERED_TOP_MASK) != LAYERED_TOP)
/** node has a neighbor below (modulo n_nodes) */
#define LAYERED_BTM_NEIGHBOR ((layered_flags & LAYERED_BTM_MASK) != LAYERED_BOTTOM)

int layered_flags = 0;
int n_layers = -1, determine_n_layers = 1;
double layer_h = 0, layer_h_i = 0;

static int btm, top;

void layered_get_mi_vector(double res[3], double a[3], double b[3])
{
  int i;

  for(i=0;i<2;i++) {
    res[i] = a[i] - b[i];
#ifdef PARTIAL_PERIODIC
    if (PERIODIC(i))
#endif
      res[i] -= dround(res[i]*box_l_i[i])*box_l[i];
  }
  res[2] = a[2] - b[2];
}

Cell *layered_position_to_cell(double pos[3])
{
  int cpos = (int)((pos[2] - my_left[2])*layer_h_i) + 1;
  if (cpos < 1) {
    if (!LAYERED_BTM_NEIGHBOR)
      cpos = 1;
    else
      return NULL;
  }
  else if (cpos > n_layers) {
    /* not periodic, but at top */
    if (!LAYERED_TOP_NEIGHBOR)
      cpos = n_layers;
    else
      return NULL;
  }
  return &(cells[cpos]);
}

void layered_topology_release()
{
  CELL_TRACE(fprintf(stderr,"%d: layered_topology_release:\n", this_node));
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
}

static void layered_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int even_odd;
  int c, n;

  if (n_nodes > 1) {
    /* more than one node => no local transfers */

    /* how many communications to do: up even/odd, down even/odd */
    n = 4;
    /* one communication missing if not periodic but on border */
    if (!LAYERED_TOP_NEIGHBOR)
      n -= 2;
    if (!LAYERED_BTM_NEIGHBOR)
      n -= 2;

    prepare_comm(comm, data_parts, n);

    /* always sending/receiving 1 cell per time step */
    for(c = 0; c < n; c++) {
      comm->comm[c].part_lists = (ParticleList**)malloc(sizeof(ParticleList *));
      comm->comm[c].n_part_lists = 1;
      comm->comm[c].mpi_comm = comm_cart;
    }

    c = 0;

    CELL_TRACE(fprintf(stderr, "%d: ghostrec new comm of size %d\n", this_node, n));
    /* downwards */
    for (even_odd = 0; even_odd < 2; even_odd++) {
      /* send */
      if (this_node % 2 == even_odd && LAYERED_BTM_NEIGHBOR) {
	comm->comm[c].type = GHOST_SEND;
	/* round 1 uses prefetched data and stores delayed data */
	if (c == 1)
	  comm->comm[c].type |= GHOST_PREFETCH | GHOST_PSTSTORE;
	comm->comm[c].node = btm;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[0];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send force to %d btmg\n", this_node, btm));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[1];

	  /* if periodic and bottom or top, send shifted */
	  comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	  if (((layered_flags & LAYERED_BTM_MASK) == LAYERED_BTM_MASK) &&
	      (data_parts & GHOSTTRANS_POSITION)) {
	    comm->data_parts |= GHOSTTRANS_POSSHFTD;
	    comm->comm[c].shift[2] = box_l[2];
	  }
	  else
	    comm->comm[c].shift[2] = 0;
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send to %d shift %f btml\n", this_node, btm, comm->comm[c].shift[2]));
	}
	c++;
      }
      /* recv. Note we test r_node as we always have to test the sender
	 as for odd n_nodes maybe we send AND receive. */
      if (top % 2 == even_odd && LAYERED_TOP_NEIGHBOR) {
	comm->comm[c].type = GHOST_RECV;
	/* round 0 prefetch send for round 1 and delay recvd data processing */
	if (c == 0)
	  comm->comm[c].type |= GHOST_PREFETCH | GHOST_PSTSTORE;
	comm->comm[c].node = top;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[n_layers];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec get force from %d topl\n", this_node, top));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[n_layers + 1];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec recv from %d topg\n", this_node, top));
	}
	c++;
      }
    }

    CELL_TRACE(fprintf(stderr, "%d: ghostrec upwards\n", this_node));
    /* upwards */
    for (even_odd = 0; even_odd < 2; even_odd++) {
      /* send */
      if (this_node % 2 == even_odd && LAYERED_TOP_NEIGHBOR) {
	comm->comm[c].type = GHOST_SEND;
	/* round 1 use prefetched data from round 0.
	   But this time there may already have been two transfers downwards */
	if (c % 2 == 1)
	  comm->comm[c].type |= GHOST_PREFETCH | GHOST_PSTSTORE;
	comm->comm[c].node = top;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[n_layers + 1];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send force to %d topg\n", this_node, top));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[n_layers];

	  /* if periodic and bottom or top, send shifted */
	  comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	  if (((layered_flags & LAYERED_TOP_MASK) == LAYERED_TOP_MASK) &&
	      (data_parts & GHOSTTRANS_POSITION)) {
	    comm->data_parts |= GHOSTTRANS_POSSHFTD;
	    comm->comm[c].shift[2] = -box_l[2];
	  }
	  else
	    comm->comm[c].shift[2] = 0;
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send to %d shift %f topl\n", this_node, top, comm->comm[c].shift[2]));
	}
	c++;
      }
      /* recv. Note we test r_node as we always have to test the sender
	 as for odd n_nodes maybe we send AND receive. */
      if (btm % 2 == even_odd && LAYERED_BTM_NEIGHBOR) {
	comm->comm[c].type = GHOST_RECV;
	/* round 0 prefetch. But this time there may already have been two transfers downwards */
	if (c % 2 == 0)
	  comm->comm[c].type |= GHOST_PREFETCH | GHOST_PSTSTORE;
	comm->comm[c].node = btm;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[1];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec get force from %d btml\n", this_node, btm));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[0];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec recv from %d btmg\n", this_node, btm));
	}
	c++;
      }
    }
  }
  else {
    /* one node => local transfers, either 2 (up and down, periodic) or zero*/

    n = (layered_flags & LAYERED_PERIODIC) ? 2 : 0;

    prepare_comm(comm, data_parts, n);

    if (n != 0) {
      /* two cells: from and to */
      for(c = 0; c < n; c++) {
	comm->comm[c].part_lists = (ParticleList**)malloc(2*sizeof(ParticleList *));
	comm->comm[c].n_part_lists = 2;
	comm->comm[c].mpi_comm = comm_cart;
	comm->comm[c].node = this_node;
      }

      c = 0;

      /* downwards */
      comm->comm[c].type = GHOST_LOCL;
      if (data_parts == GHOSTTRANS_FORCE) {
	comm->comm[c].part_lists[0] = &cells[0];
	comm->comm[c].part_lists[1] = &cells[n_layers];
      }
      else {
	comm->comm[c].part_lists[0] = &cells[1];
	comm->comm[c].part_lists[1] = &cells[n_layers + 1];
	/* here it is periodic */
	if (data_parts & GHOSTTRANS_POSITION)
	  comm->data_parts |= GHOSTTRANS_POSSHFTD;
	comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	comm->comm[c].shift[2] = box_l[2];
      }
      c++;

      /* upwards */
      comm->comm[c].type = GHOST_LOCL;
      if (data_parts == GHOSTTRANS_FORCE) {
	comm->comm[c].part_lists[0] = &cells[n_layers + 1];
	comm->comm[c].part_lists[1] = &cells[1];
      }
      else {
	comm->comm[c].part_lists[0] = &cells[n_layers];
	comm->comm[c].part_lists[1] = &cells[0];
	/* here it is periodic */
	if (data_parts & GHOSTTRANS_POSITION)
	  comm->data_parts |= GHOSTTRANS_POSSHFTD;
	comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	comm->comm[c].shift[2] = -box_l[2];
      }
    }
  }
}

void layered_topology_init(CellPList *old)
{
  Particle *part;
  int c, p, np;

  CELL_TRACE(fprintf(stderr, "%d: layered_topology_init, %d old particle lists\n", this_node, old->n));

  cell_structure.type = CELL_STRUCTURE_LAYERED;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = layered_position_to_cell;

  /* check node grid. All we can do is 1x1xn. */
  if (node_grid[0] != 1 || node_grid[1] != 1) {
    char *errtxt = runtime_error(128 + ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{016 selected node grid is not suitable for layered cell structure (needs 1x1x%d grid)} ", n_nodes);
    node_grid[0] = node_grid[1] = 1;
    node_grid[2] = n_nodes;
  }

  if (this_node == 0 && determine_n_layers) {
    if (max_range > 0) {
      n_layers = (int)floor(local_box_l[2]/max_range);
      if (n_layers < 1) {
	char *errtxt = runtime_error(128 + 2*ES_DOUBLE_SPACE);
	ERROR_SPRINTF(errtxt, "{017 layered: maximal interaction range %g larger than local box length %g} ", max_range, local_box_l[2]);
	n_layers = 1;
      }
      if (n_layers > max_num_cells - 2)
	n_layers = imax(max_num_cells - 2, 0);
    }
    else
      n_layers = 1;
  }
  MPI_Bcast(&n_layers, 1, MPI_INT, 0, comm_cart);

  top = this_node + 1;
  if (top == n_nodes && (layered_flags & LAYERED_PERIODIC))
    top = 0;
  btm = this_node - 1;
  if (btm == -1 && (layered_flags & LAYERED_PERIODIC))
    btm = n_nodes - 1;

  layer_h = local_box_l[2]/(double)(n_layers);
  layer_h_i = 1/layer_h;

  if (layer_h < max_range) {
    char *errtxt = runtime_error(128 + 2*ES_DOUBLE_SPACE);
    ERROR_SPRINTF(errtxt, "{018 layered: maximal interaction range %g larger than layer height %g} ", max_range, layer_h);
  }

  /* check wether node is top and/or bottom */
  layered_flags = 0;
  if (this_node == 0)
    layered_flags |= LAYERED_BOTTOM;
  if (this_node == n_nodes - 1)
    layered_flags |= LAYERED_TOP;

  if (PERIODIC(2))
    layered_flags |= LAYERED_PERIODIC;

  CELL_TRACE(fprintf(stderr, "%d: layered_flags tn %d bn %d \n", this_node, LAYERED_TOP_NEIGHBOR, LAYERED_BTM_NEIGHBOR));

  /* allocate cells and mark them */
  realloc_cells(n_layers + 2);
  realloc_cellplist(&local_cells, local_cells.n = n_layers);
  for (c = 0; c < n_layers; c++)
    local_cells.cell[c] = &cells[c + 1];
  realloc_cellplist(&ghost_cells, ghost_cells.n = 2);
  ghost_cells.cell[0] = &cells[0];
  ghost_cells.cell[1] = &cells[n_layers + 1];

  /* create communicators */
  layered_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);
  layered_prepare_comm(&cell_structure.exchange_ghosts_comm,     GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  layered_prepare_comm(&cell_structure.update_ghost_pos_comm,    GHOSTTRANS_POSITION);
  layered_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++) {
      Cell *nc = layered_position_to_cell(part[p].r.p);
      /* particle does not belong to this node. Just stow away
	 somewhere for the moment */
      if (nc == NULL)
	nc = local_cells.cell[0];
      append_unindexed_particle(nc, &part[p]);
    }
  }
  for (c = 1; c <= n_layers; c++)
    update_local_particles(&cells[c]);

  CELL_TRACE(fprintf(stderr, "%d: layered_topology_init done\n", this_node));
}
 
static void layered_append_particles(ParticleList *pl, ParticleList *up, ParticleList *dn)
{
  int p;

  CELL_TRACE(fprintf(stderr, "%d: sorting in %d\n", this_node, pl->n));
  for(p = 0; p < pl->n; p++) {
    fold_position(pl->part[p].r.p, pl->part[p].l.i);
    if (LAYERED_BTM_NEIGHBOR && pl->part[p].r.p[2] < my_left[2]) {
      CELL_TRACE(fprintf(stderr, "%d: leaving part %d for node below\n", this_node, pl->part[p].p.identity));
      move_indexed_particle(dn, pl, p);
    }
    else if (LAYERED_TOP_NEIGHBOR && pl->part[p].r.p[2] >= my_right[2]) {
      CELL_TRACE(fprintf(stderr, "%d: leaving part %d for node above\n", this_node, pl->part[p].p.identity));
      move_indexed_particle(up, pl, p);
    }
    else
      move_indexed_particle(layered_position_to_cell(pl->part[p].r.p), pl, p);
    /* same particle again, as this is now a new one */
    if (p < pl->n) p--;
  }
  CELL_TRACE(fprintf(stderr, "%d: left over %d\n", this_node, pl->n));
}

void layered_exchange_and_sort_particles(int global_flag)
{
  Particle *part;
  Cell *nc, *oc;
  int c, p, flag, redo;
  ParticleList send_buf_dn, send_buf_up;
  ParticleList recv_buf;

  CELL_TRACE(fprintf(stderr, "%d:layered exchange and sort %d\n", this_node, global_flag));

  init_particlelist(&send_buf_dn);
  init_particlelist(&send_buf_up);

  init_particlelist(&recv_buf);

  /* sort local particles */
  for (c = 1; c <= n_layers; c++) {
    oc = &cells[c];

    for (p = 0; p < oc->n; p++) {
      part = &oc->part[p];

      if (n_nodes != 1 && LAYERED_BTM_NEIGHBOR && part->r.p[2] < my_left[2]) {
	CELL_TRACE(fprintf(stderr, "%d: send part %d down\n", this_node, part->p.identity));
	move_indexed_particle(&send_buf_dn, oc, p);
	if (p < oc->n) p--;
      }
      else if (n_nodes != 1 && LAYERED_TOP_NEIGHBOR && part->r.p[2] >= my_right[2]) {
	CELL_TRACE(fprintf(stderr, "%d: send part %d up\n", this_node, part->p.identity));
	move_indexed_particle(&send_buf_up, oc, p);
	if (p < oc->n) p--;
      }
      else {
	/* particle stays here. Fold anyways to get x,y correct */
	fold_position(part->r.p, part->l.i);
	nc = layered_position_to_cell(part->r.p);
	if (nc != oc) {
	  move_indexed_particle(nc, oc, p);
	  if (p < oc->n) p--;
	}
      }
    }
  }

  for (;;) {
    /* transfer */
    if (n_nodes > 1) {
      if (this_node % 2 == 0) {
	/* send down */
	if (LAYERED_BTM_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: send dn\n",this_node));
	  send_particles(&send_buf_dn, btm);
	}
	if (LAYERED_TOP_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: recv up\n",this_node));
	  recv_particles(&recv_buf, top);
	}
	/* send up */
	if (LAYERED_TOP_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: send up\n",this_node));
	  send_particles(&send_buf_up, top);
	}
	if (LAYERED_BTM_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: recv dn\n",this_node));
	  recv_particles(&recv_buf, btm);
	}
      }
      else {
	if (LAYERED_TOP_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: recv up\n",this_node));
	  recv_particles(&recv_buf, top);
	}
	if (LAYERED_BTM_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: send dn\n",this_node));
	  send_particles(&send_buf_dn, btm);
	}
	if (LAYERED_BTM_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: recv dn\n",this_node));
	  recv_particles(&recv_buf, btm);
	}
	if (LAYERED_TOP_NEIGHBOR) {
	  CELL_TRACE(fprintf(stderr,"%d: send up\n",this_node));
	  send_particles(&send_buf_up, top);
	}
      }
    }
    else {
      if (recv_buf.n != 0 || send_buf_dn.n != 0 || send_buf_up.n != 0) {
	fprintf(stderr, "1 node but transfer buffers are not empty. send up %d, down %d, recv %d\n",
		send_buf_up.n, send_buf_dn.n, recv_buf.n);
	errexit();
      }
    }
    layered_append_particles(&recv_buf, &send_buf_up, &send_buf_dn);

    /* handshake redo */
    flag = (send_buf_up.n != 0 || send_buf_dn.n != 0);
    
    CELL_TRACE(if (flag) fprintf(stderr, "%d: requesting another exchange round\n", this_node));

    if (global_flag == CELL_GLOBAL_EXCHANGE) {
      MPI_Allreduce(&flag, &redo, 1, MPI_INT, MPI_MAX, comm_cart);
      if (!redo)
	break;
      CELL_TRACE(fprintf(stderr, "%d: another exchange round\n", this_node));
    }
    else {
      if (flag) {
	char *errtxt = runtime_error(128 + ES_DOUBLE_SPACE);
	ERROR_SPRINTF(errtxt,"{019 layered_exchange_and_sort_particles: particle moved more than one cell} ");

	/* sort left over particles into border cells */
	CELL_TRACE(fprintf(stderr, "%d: emergency sort\n", this_node));
	while (send_buf_up.n > 0)
	  move_indexed_particle(&cells[1], &send_buf_up, 0);
	while (send_buf_dn.n > 0)
	  move_indexed_particle(&cells[n_layers], &send_buf_dn, 0);
      }
      break;
    }
  }

  realloc_particlelist(&recv_buf, 0);
}

/** nonbonded and bonded force calculation using the verlet list */
void layered_calculate_ia()
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];
 
  CELL_TRACE(fprintf(stderr, "%d: rebuild_v=%d\n", this_node, rebuild_verletlist));

  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    cellb = &cells[c-1];
    pb    = cellb->part;
    npb   = cellb->n;

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      if (rebuild_verletlist)
	memcpy(p1->l.p_old, p1->r.p, 3*sizeof(double));

      add_bonded_force(p1);
#ifdef CONSTRAINTS
      add_constraints_forces(p1);
#endif
      add_external_potential_forces(p1);

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
	layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
	dist2 = sqrlen(d);
#ifdef EXCLUSIONS
	if (do_nonbonded(p1, &pl[j]))
#endif
	  add_non_bonded_pair_force(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
	layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
	dist2 = sqrlen(d);
#ifdef EXCLUSIONS
	if (do_nonbonded(p1, &pb[j]))
#endif
	  add_non_bonded_pair_force(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
  rebuild_verletlist = 0;
}

void layered_calculate_energies()
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];
 
  CELL_TRACE(fprintf(stderr, "%d: rebuild_v=%d\n", this_node, rebuild_verletlist));

  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    cellb = &cells[c-1];
    pb    = cellb->part;
    npb   = cellb->n;

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      if (rebuild_verletlist)
	memcpy(p1->l.p_old, p1->r.p, 3*sizeof(double));

      add_kinetic_energy(p1);

      add_bonded_energy(p1);
#ifdef CONSTRAINTS
      add_constraints_energy(p1);
#endif
      add_external_potential_energy(p1);

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
	layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
	dist2 = sqrlen(d);
#ifdef EXCLUSIONS
	if (do_nonbonded(p1, &pl[j]))
#endif
	  add_non_bonded_pair_energy(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
	layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
	dist2 = sqrlen(d);
#ifdef EXCLUSIONS
	if (do_nonbonded(p1, &pb[j]))
#endif
	  add_non_bonded_pair_energy(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
  rebuild_verletlist = 0;
}

void layered_calculate_virials(int v_comp)
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];
 
  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    cellb = &cells[c-1];
    pb    = cellb->part;
    npb   = cellb->n;

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      if (rebuild_verletlist)
	memcpy(p1->l.p_old, p1->r.p, 3*sizeof(double));

      add_kinetic_virials(p1,v_comp);

      add_bonded_virials(p1);
#ifdef BOND_ANGLE_OLD
      add_three_body_bonded_stress(p1);
#endif
#ifdef BOND_ANGLE
      add_three_body_bonded_stress(p1);
#endif

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
	layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
	dist2 = sqrlen(d);
#ifdef EXCLUSIONS
	if (do_nonbonded(p1, &pl[j]))
#endif
	  add_non_bonded_pair_virials(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
	layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
	dist2 = sqrlen(d);
#ifdef EXCLUSIONS
	if (do_nonbonded(p1, &pb[j]))
#endif
	  add_non_bonded_pair_virials(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
  rebuild_verletlist = 0;
}
