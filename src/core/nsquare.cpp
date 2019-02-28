/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
 *
 *  Implementation of nsquare.hpp.
 */

#include "nsquare.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "debug.hpp"
#include "ghosts.hpp"
#include "utils.hpp"

#include <cstring>
#include <mpi.h>

Cell *local;

Cell *nsq_position_to_cell(const Vector3d &pos) { return local; }
int nsq_position_to_node(const Vector3d &) { return this_node; }

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

void nsq_topology_init() {
  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.position_to_node = nsq_position_to_node;
  cell_structure.position_to_cell = nsq_position_to_cell;

  realloc_cells(n_nodes);

  /* mark cells */
  local = &cells[this_node];
  realloc_cellplist(&local_cells, local_cells.n = 1);
  local_cells.cell[0] = local;

  realloc_cellplist(&ghost_cells, ghost_cells.n = n_nodes - 1);
  int c = 0;
  for (int n = 0; n < n_nodes; n++)
    if (n != this_node)
      ghost_cells.cell[c++] = &cells[n];

  std::vector<Cell *> red_neighbors;
  std::vector<Cell *> black_neighbors;

  /* distribute force calculation work  */
  for (int n = 0; n < n_nodes; n++) {
    auto const diff = n - this_node;
    /* simple load balancing formula. Basically diff % n, where n >= n_nodes, n
       odd.
       The node itself is also left out, as it is treated differently */
    if (diff == 0) {
      continue;
    }

    if (((diff > 0 && (diff % 2) == 0) || (diff < 0 && ((-diff) % 2) == 1))) {
      CELL_TRACE(
          fprintf(stderr, "%d: doing interactions with %d\n", this_node, n));
      red_neighbors.push_back(&cells.at(n));
    } else {
      black_neighbors.push_back(&cells.at(n));
    }
  }

  local->m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);

  /* create communicators */
  nsq_prepare_comm(&cell_structure.ghost_cells_comm, GHOSTTRANS_PARTNUM);
  nsq_prepare_comm(&cell_structure.exchange_ghosts_comm,
                   GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  nsq_prepare_comm(&cell_structure.update_ghost_pos_comm, GHOSTTRANS_POSITION);
  nsq_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

  /* here we just decide what to transfer where */
  if (n_nodes > 1) {
    for (int n = 0; n < n_nodes; n++) {
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
}

ParticleList nsq_balance_particles(int global_flag) {
    /* we don't have the concept of neighbors, and therefore don't need that.
     However, if global particle changes happen, we might want to rebalance. */
  if (global_flag != CELL_GLOBAL_EXCHANGE)
    return {};

  std::vector<int> ppnode;
  boost::mpi::all_gather(comm_cart, cells_get_n_particles(), ppnode);

  /* minimal difference between node shares */
  int minshare = n_part / n_nodes;
  int maxshare = minshare + 1;

  CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: load %d-%d\n",
                     this_node, minshare, maxshare));

  ParticleList new_particles;
  std::vector<int> removed_particles;
  for (;;) {
    /* find node with most excessive particles */
    auto const minmax = std::minmax_element(ppnode.begin(), ppnode.end());
    /* Nothing to do */
    if(minmax.first == minmax.second) {
      break;
    }

    auto const surplus = *minmax.second - minshare;
    auto const s_node = std::distance(ppnode.begin(), minmax.second);
    auto const lack = maxshare - *minmax.first;
    auto const l_node = std::distance(ppnode.begin(), minmax.first);

    CELL_TRACE(fprintf(stderr,
                       "%d: nsq_balance_particles: lack %d on node %d\n",
                       this_node, lack, l_node));

    /* exit if all nodes load is withing min and max share */
    if (lack <= 1 && surplus <= 1)
      break;

    auto const transfer = lack < surplus ? lack : surplus;

    if (s_node == this_node) {
      ParticleList send_buf;

      realloc_particlelist(&send_buf, send_buf.n = transfer);

      {
        int i = 0;
        for (; i < std::min(transfer, new_particles.n); i++) {
          send_buf.part[i] = std::move(new_particles.part[--new_particles.n]);
        }
        for (; i < transfer; i++) {
          send_buf.part[i] = std::move(local->part[--local->n]);
        }
      }

      realloc_particlelist(local, local->n);

      send_particles(&send_buf, l_node);
    } else if (l_node == this_node) {
      ParticleList recv_buf{};

      recv_particles(&recv_buf, s_node);
      for (int i = 0; i < recv_buf.n; i++) {
        append_particle(&new_particles, std::move(recv_buf.part[i]));
      }

      realloc_particlelist(&recv_buf, 0);
    }
    ppnode[s_node] -= transfer;
    ppnode[l_node] += transfer;
  }
  CELL_TRACE(fprintf(stderr, "%d: nsq_balance_particles: done\n", this_node));

  for (int i = 0; i < new_particles.n; i++) {
    append_particle(local, std::move(new_particles.part[i]));
  }

  return {};
}
