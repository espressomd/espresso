/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *
 *  Implementation of nsquare.hpp.
 */

#include "nsquare.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "ghosts.hpp"

#include <mpi.h>

Cell *local;

static Cell *nsq_id_to_cell(int id) {
  return ((id % n_nodes) == this_node) ? local : nullptr;
}

void nsq_topology_release() {
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
  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.particle_to_cell = [](const Particle &p) {
    return nsq_id_to_cell(p.identity());
  };

  /* This cell system supports any range that is compatible with
   * the geometry. */
  for (int i = 0; i < 3; i++) {
    cell_structure.max_range[i] = box_geo.periodic(i)
                                      ? 0.5 * box_geo.length()[i]
                                      : std::numeric_limits<double>::infinity();
  }

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

  /* copy particles */
  for (int c = 0; c < old->n; c++) {
    auto part = old->cell[c]->part;
    auto np = old->cell[c]->n;
    for (int p = 0; p < np; p++)
      append_unindexed_particle(local, std::move(part[p]));
  }
  update_local_particles(local);
}

void nsq_exchange_particles(int global_flag, ParticleList *displaced_parts) {
  if (not global_flag) {
    assert(displaced_parts->n == 0);
    return;
  }

  /* Sort displaced particles by the node they belong to. */
  std::vector<std::vector<Particle>> send_buf(n_nodes);
  for (auto &p : Utils::make_span(displaced_parts->part, displaced_parts->n)) {
    auto const target_node = (p.identity() % n_nodes);
    send_buf.at(target_node).emplace_back(std::move(p));
  }
  realloc_particlelist(displaced_parts, displaced_parts->n = 0);

  /* Exchange particles */
  std::vector<std::vector<Particle>> recv_buf(n_nodes);
  boost::mpi::all_to_all(comm_cart, send_buf, recv_buf);

  /* Add new particles belonging to this node */
  for (auto &parts : recv_buf) {
    for (auto &p : parts) {
      append_indexed_particle(local, std::move(p));
    }
  }
}
