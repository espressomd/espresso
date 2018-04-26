/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file nsquare.cpp
 *
 *  Implementation of  \ref nsquare.hpp "nsquare.h".
 */

#include "nsquare.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "ghosts.hpp"
#include "utils.hpp"
#include "debug.hpp"

#include <cstring>
#include <mpi.h>

Cell *local;

Cell *nsq_position_to_cell(double pos[3]) { return local; }

void nsq_topology_release() {
  CELL_TRACE(fprintf(stderr, "%d: nsq_topology_release:\n", this_node));
  /* free ghost cell pointer list */
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
}

static void nsq_prepare_comm(GhostCommunicator *comm, int data_parts) {
  int n;
  /* no need for comm for only 1 node */
  if (n_nodes == 1) {
    prepare_comm(comm, data_parts, 0);
    return;
  }

  prepare_comm(comm, data_parts, n_nodes);
  /* every node has its dedicated comm step */
  for (n = 0; n < n_nodes; n++) {
    comm->comm[n].part_lists = (Cell **)Utils::malloc(sizeof(Cell *));
    comm->comm[n].part_lists[0] = &cells[n];
    comm->comm[n].n_part_lists = 1;
    comm->comm[n].node = n;
    comm->comm[n].mpi_comm = comm_cart;
  }
}

void nsq_topology_init(CellPList *old) {
  Particle *part;
  int n, p, np, diff;

  CELL_TRACE(fprintf(stderr, "%d: nsq_topology_init, %d\n", this_node, old->n));

  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = nsq_position_to_cell;

  realloc_cells(n_nodes);

  /* mark cells */
  local = &cells[this_node];
  realloc_cellplist(&local_cells, local_cells.n = 1);
  local_cells.cell[0] = local;

  cells[this_node].m_neighbors.clear();

  realloc_cellplist(&ghost_cells, ghost_cells.n = n_nodes - 1);
  auto c = 0;
  for (n = 0; n < n_nodes; n++)
    if (n != this_node)
      ghost_cells.cell[c++] = &cells[n];

  /* distribute force calculation work  */

  for (n = 0; n < n_nodes; n++) {
    diff = n - this_node;
    /* simple load balancing formula. Basically diff % n, where n >= n_nodes, n
       odd.
       The node itself is also left out, as it is treated differently */
    if (((diff > 0 && (diff % 2) == 0) || (diff < 0 && ((-diff) % 2) == 1))) {
      CELL_TRACE(
          fprintf(stderr, "%d: doing interactions with %d\n", this_node, n));
      cells[this_node].m_neighbors.push_back(std::ref(cells[n]));
    }
  }

  /* create communicators */
  nsq_prepare_comm(&cell_structure.ghost_cells_comm, GHOSTTRANS_PARTNUM);
  nsq_prepare_comm(&cell_structure.exchange_ghosts_comm,
                   GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  nsq_prepare_comm(&cell_structure.update_ghost_pos_comm, GHOSTTRANS_POSITION);
  nsq_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

  /* here we just decide what to transfer where */
  if (n_nodes > 1) {
    for (n = 0; n < n_nodes; n++) {
      /* use the prefetched send buffers. Node 0 transmits first and never
       * prefetches. */
      if (this_node == 0 || this_node != n) {
        cell_structure.ghost_cells_comm.comm[n].type = GHOST_BCST;
        cell_structure.exchange_ghosts_comm.comm[n].type = GHOST_BCST;
        cell_structure.update_ghost_pos_comm.comm[n].type = GHOST_BCST;
      } else {
        cell_structure.ghost_cells_comm.comm[n].type =
            GHOST_BCST | GHOST_PREFETCH;
        cell_structure.exchange_ghosts_comm.comm[n].type =
            GHOST_BCST | GHOST_PREFETCH;
        cell_structure.update_ghost_pos_comm.comm[n].type =
            GHOST_BCST | GHOST_PREFETCH;
      }
      cell_structure.collect_ghost_force_comm.comm[n].type = GHOST_RDCE;
    }
    /* first round: all nodes except the first one prefetch their send data */
    if (this_node != 0) {
      cell_structure.ghost_cells_comm.comm[0].type |= GHOST_PREFETCH;
      cell_structure.exchange_ghosts_comm.comm[0].type |= GHOST_PREFETCH;
      cell_structure.update_ghost_pos_comm.comm[0].type |= GHOST_PREFETCH;
    }
  }

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np = old->cell[c]->n;
    for (p = 0; p < np; p++)
      append_unindexed_particle(local, std::move(part[p]));
  }
  update_local_particles(local);
}

void nsq_balance_particles(int global_flag) {
  int i, n, surplus, s_node, tmp, lack, l_node, transfer;

  /* Refold positions in any case, this is always safe. */
  for (auto &p : local_cells.particles()) {
    fold_position(p.r.p, p.l.i);
  }

  /* we don't have the concept of neighbors, and therefore don't need that.
     However, if global particle changes happen, we might want to rebalance. */
  if (global_flag != CELL_GLOBAL_EXCHANGE)
    return;

  int pp = cells_get_n_particles();
  int *ppnode = (int *)Utils::malloc(n_nodes * sizeof(int));
  /* minimal difference between node shares */
  int minshare = n_part / n_nodes;
  int maxshare = minshare + 1;

  CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: load %d-%d\n",
                     this_node, minshare, maxshare));

  MPI_Allgather(&pp, 1, MPI_INT, ppnode, 1, MPI_INT, comm_cart);
  for (;;) {
    /* find node with most excessive particles */
    surplus = -1;
    s_node = -1;
    for (n = 0; n < n_nodes; n++) {
      tmp = ppnode[n] - minshare;
      CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: node %d has %d\n",
                         this_node, n, ppnode[n]));
      if (tmp > surplus) {
        surplus = tmp;
        s_node = n;
      }
    }
    CELL_TRACE(fprintf(stderr,
                       "%d: nsq_balance_particles: excess %d on node %d\n",
                       this_node, surplus, s_node));

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
    CELL_TRACE(fprintf(stderr,
                       "%d: nsq_balance_particles: lack %d on node %d\n",
                       this_node, lack, l_node));

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
        send_buf.part[i] = std::move(local->part[--local->n]);
      }
      realloc_particlelist(local, local->n);
      update_local_particles(local);

      send_particles(&send_buf, l_node);
#ifdef ADDITIONAL_CHECKS
      check_particle_consistency();
#endif
    }
    else if (l_node == this_node) {
      ParticleList recv_buf{};

      recv_particles(&recv_buf, s_node);
      for (int i = 0; i < recv_buf.n; i++) {
        append_indexed_particle(local, std::move(recv_buf.part[i]));
      }

      realloc_particlelist(&recv_buf, 0);

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
