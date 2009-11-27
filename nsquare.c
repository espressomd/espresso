// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file nsquare.c
 *
 *  Implementation of  \ref nsquare.h "nsquare.h".
 */

#include <mpi.h>
#include <string.h>
#include "utils.h"
#include "nsquare.h"
#include "communication.h"
#include "ghosts.h"
#include "forces.h"
#include "pressure.h"
#include "energy.h"

Cell *local;
CellPList me_do_ghosts;

Cell *nsq_position_to_cell(double pos[3])
{
  return local;
}

void nsq_topology_release()
{
  CELL_TRACE(fprintf(stderr,"%d: nsq_topology_release:\n",this_node));
  /* free ghost cell pointer list */
  realloc_cellplist(&me_do_ghosts, 0);
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
}

static void nsq_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int n;
  /* no need for comm for only 1 node */
  if (n_nodes == 1) {
    prepare_comm(comm, data_parts, 0);
    return;
  }

  prepare_comm(comm, data_parts, n_nodes);
  /* every node has its dedicated comm step */
  for(n = 0; n < n_nodes; n++) {
    comm->comm[n].part_lists = malloc(sizeof(ParticleList *));
    comm->comm[n].part_lists[0] = &cells[n];
    comm->comm[n].n_part_lists = 1;
    comm->comm[n].node = n;
    comm->comm[n].mpi_comm = MPI_COMM_WORLD;
  }
}

void nsq_topology_init(CellPList *old)
{
  Particle *part;
  int n, c, p, np, ntodo, diff;

  CELL_TRACE(fprintf(stderr, "%d: nsq_topology_init, %d\n", this_node, old->n));

  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = nsq_position_to_cell;

  realloc_cells(n_nodes);

  /* mark cells */
  local = &cells[this_node];
  realloc_cellplist(&local_cells, local_cells.n = 1);
  local_cells.cell[0] = local;

  realloc_cellplist(&ghost_cells, ghost_cells.n = n_nodes - 1);
  c = 0;
  for (n = 0; n < n_nodes; n++)
    if (n != this_node)
      ghost_cells.cell[c++] = &cells[n];

  /* distribute force calculation work  */
  ntodo = (n_nodes + 3)/2;
  init_cellplist(&me_do_ghosts);
  realloc_cellplist(&me_do_ghosts, ntodo);
  for (n = 0; n < n_nodes; n++) {
    diff = n - this_node;
    /* simple load balancing formula. Basically diff % n, where n >= n_nodes, n odd.
       The node itself is also left out, as it is treated differently */
    if (((diff > 0 && (diff % 2) == 0) ||
	 (diff < 0 && ((-diff) % 2) == 1))) {
      CELL_TRACE(fprintf(stderr, "%d: doing interactions with %d\n", this_node, n));
      me_do_ghosts.cell[me_do_ghosts.n++] = &cells[n];
    }
  }

  /* create communicators */
  nsq_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);
  nsq_prepare_comm(&cell_structure.exchange_ghosts_comm,     GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  nsq_prepare_comm(&cell_structure.update_ghost_pos_comm,    GHOSTTRANS_POSITION);
  nsq_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

  /* here we just decide what to transfer where */
  if (n_nodes > 1) {
    for (n = 0; n < n_nodes; n++) {
      /* use the prefetched send buffers. Node 0 transmits first and never prefetches. */
      if (this_node == 0 || this_node != n) {
	cell_structure.ghost_cells_comm.comm[n].type         = GHOST_BCST;
	cell_structure.exchange_ghosts_comm.comm[n].type     = GHOST_BCST;
	cell_structure.update_ghost_pos_comm.comm[n].type    = GHOST_BCST;
      }
      else {
	cell_structure.ghost_cells_comm.comm[n].type         = GHOST_BCST | GHOST_PREFETCH;
	cell_structure.exchange_ghosts_comm.comm[n].type     = GHOST_BCST | GHOST_PREFETCH;
	cell_structure.update_ghost_pos_comm.comm[n].type    = GHOST_BCST | GHOST_PREFETCH;
      }
      cell_structure.collect_ghost_force_comm.comm[n].type = GHOST_RDCE;
    }
    /* first round: all nodes except the first one prefetch their send data */
    if (this_node != 0) {
      cell_structure.ghost_cells_comm.comm[0].type         |= GHOST_PREFETCH;
      cell_structure.exchange_ghosts_comm.comm[0].type     |= GHOST_PREFETCH;
      cell_structure.update_ghost_pos_comm.comm[0].type    |= GHOST_PREFETCH;
    }
  }

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++)
      append_unindexed_particle(local, &part[p]);
  }
  update_local_particles(local);
}

void nsq_balance_particles()
{
  int i, n, surplus, s_node, tmp, lack, l_node, transfer;

  int pp = cells_get_n_particles();
  int *ppnode = malloc(n_nodes*sizeof(int));
  /* minimal difference between node shares */
  int minshare = n_total_particles/n_nodes;
  int maxshare = minshare + 1;

  CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: load %d-%d\n", this_node, minshare, maxshare));

  MPI_Allgather(&pp, 1, MPI_INT, ppnode, 1, MPI_INT, MPI_COMM_WORLD);
  for (;;) {
    /* find node with most excessive particles */
    surplus = -1;
    s_node = -1;
    for (n = 0; n < n_nodes; n++) {
      tmp = ppnode[n] - minshare;
      CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: node %d has %d\n", this_node, n, ppnode[n]));
      if (tmp > surplus) {
	surplus = tmp;
	s_node = n;
      }
    }
    CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: excess %d on node %d\n", this_node, surplus, s_node));

    /* find node with most lacking particles */
    lack = -1;
    l_node = -1;
    for (n = 0; n < n_nodes; n++) {
      tmp = maxshare - ppnode[n];
      if (tmp > lack) {
	lack = tmp;
	l_node = n;
      }
    }
    CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: lack %d on node %d\n", this_node, lack, l_node));

    /* should not happen: minshare or maxshare wrong or more likely,
       the algorithm */
    if (s_node == -1 || l_node == -1) {
      fprintf(stderr, "%d: Particle load balancing failed\n", this_node);
      break;
    }

    /* exit if all nodes load is withing min and max share */
    if (lack <= 1 && surplus <= 1)
      break;

    transfer = lack < surplus ? lack : surplus;

    if (s_node == this_node) {
      ParticleList send_buf;
      init_particlelist(&send_buf);
      realloc_particlelist(&send_buf, send_buf.n = transfer);
      for (i = 0; i < transfer; i++) {
	memcpy(&send_buf.part[i], &local->part[--local->n], sizeof(Particle));
      }
      realloc_particlelist(local, local->n);
      update_local_particles(local);

      send_particles(&send_buf, l_node);
#ifdef ADDITIONAL_CHECKS
      check_particle_consistency();
#endif
    }
    else if (l_node == this_node) {
      recv_particles(local, s_node);
#ifdef ADDITIONAL_CHECKS
      check_particle_consistency();
#endif
    }
    ppnode[s_node] -= transfer;
    ppnode[l_node] += transfer;
  }
  CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: done\n", this_node));

  free(ppnode);
}

/** nonbonded and bonded force calculation using the verlet list */
void nsq_calculate_ia()
{
  Particle *partl, *partg;
  Particle *pt1, *pt2;
  int p, p2, npl, npg, c;
  double d[3], dist2, dist;

  npl   = local->n;
  partl = local->part;

  /* calculate bonded interactions and non bonded node-node */
  for (p = 0; p < npl; p++) {
    pt1 = &partl[p];
    add_bonded_force(pt1);
#ifdef CONSTRAINTS
    add_constraints_forces(pt1);
#endif

    /* other particles, same node */
    for (p2 = p + 1; p2 < npl; p2++) {
      pt2 = &partl[p2];
      get_mi_vector(d, pt1->r.p, pt2->r.p);
      dist2 = sqrlen(d);
      dist = sqrt(dist2);
#ifdef EXCLUSIONS
      if (do_nonbonded(pt1, pt2))
#endif
	add_non_bonded_pair_force(pt1, pt2, d, dist, dist2);
    }

    /* calculate with my ghosts */
    for (c = 0; c < me_do_ghosts.n; c++) {
      npg   = me_do_ghosts.cell[c]->n;
      partg = me_do_ghosts.cell[c]->part;

      for (p2 = 0; p2 < npg; p2++) {
	pt2 = &partg[p2];
	get_mi_vector(d, pt1->r.p, pt2->r.p);
	dist2 = sqrlen(d);
	dist = sqrt(dist2);
#ifdef EXCLUSIONS
	if (do_nonbonded(pt1, pt2))
#endif
	  add_non_bonded_pair_force(pt1, pt2, d, dist, dist2);
      }
    }
  }
}

void nsq_calculate_energies()
{
  Particle *partl, *partg;
  Particle *pt1, *pt2;
  int p, p2, npl, npg, c;
  double d[3], dist2, dist;

  npl   = local->n;
  partl = local->part;

  /* calculate bonded interactions and non bonded node-node */
  for (p = 0; p < npl; p++) {
    pt1 = &partl[p];
    add_kinetic_energy(pt1);
    add_bonded_energy(pt1);
#ifdef CONSTRAINTS
    add_constraints_energy(pt1);
#endif

    /* other particles, same node */
    for (p2 = p + 1; p2 < npl; p2++) {
      pt2 = &partl[p2];
      get_mi_vector(d, pt1->r.p, pt2->r.p);
      dist2 = sqrlen(d);
      dist = sqrt(dist2);
#ifdef EXCLUSIONS
      if (do_nonbonded(pt1, pt2))
#endif
	add_non_bonded_pair_energy(pt1, pt2, d, dist, dist2);
    }

    /* calculate with my ghosts */
    for (c = 0; c < me_do_ghosts.n; c++) {
      npg   = me_do_ghosts.cell[c]->n;
      partg = me_do_ghosts.cell[c]->part;

      for (p2 = 0; p2 < npg; p2++) {
	pt2 = &partg[p2];
	get_mi_vector(d, pt1->r.p, pt2->r.p);
	dist2 = sqrlen(d);
	dist = sqrt(dist2);
#ifdef EXCLUSIONS
	if (do_nonbonded(pt1, pt2))
#endif
	  add_non_bonded_pair_energy(pt1, pt2, d, dist, dist2);
      }
    }
  }
}

void nsq_calculate_virials()
{
  Particle *partl, *partg;
  Particle *pt1, *pt2;
  int p, p2, npl, npg, c;
  double d[3], dist2, dist;

  npl   = local->n;
  partl = local->part;

  /* calculate bonded interactions and non bonded node-node */
  for (p = 0; p < npl; p++) {
    pt1 = &partl[p];
    add_kinetic_virials(pt1,0);
    add_bonded_virials(pt1);
#ifdef BOND_ANGLE_COSINE
    add_three_body_bonded_stress(pt1);
#endif

    /* other particles, same node */
    for (p2 = p + 1; p2 < npl; p2++) {
      pt2 = &partl[p2];
      get_mi_vector(d, pt1->r.p, pt2->r.p);
      dist2 = sqrlen(d);
      dist = sqrt(dist2);
#ifdef EXCLUSIONS
      if (do_nonbonded(pt1, pt2))
#endif
	add_non_bonded_pair_virials(pt1, pt2, d, dist, dist2);
    }

    /* calculate with my ghosts */
    for (c = 0; c < me_do_ghosts.n; c++) {
      npg   = me_do_ghosts.cell[c]->n;
      partg = me_do_ghosts.cell[c]->part;

      for (p2 = 0; p2 < npg; p2++) {
	pt2 = &partg[p2];
	get_mi_vector(d, pt1->r.p, pt2->r.p);
	dist2 = sqrlen(d);
	dist = sqrt(dist2);
#ifdef EXCLUSIONS
	if (do_nonbonded(pt1, pt2))
#endif
	  add_non_bonded_pair_virials(pt1, pt2, d, dist, dist2);
      }
    }
  }
}
