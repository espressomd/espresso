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

static std::vector<GhostCommunication>
nsq_prepare_comm(boost::mpi::communicator const &comm) {
  std::vector<GhostCommunication> comms;

  /* every node has its dedicated comm step */
  for (int n = 0; n < comm.size(); n++) {
    GhostCommunication c;
    c.type = GHOST_BCST;

    if (comm.rank() == n) {
      c.send_lists.push_back(&cells[n]);
    } else {
      c.recv_lists.push_back(&cells[n]);
    }

    c.send_to = n;

    comms.emplace_back(std::move(c));
  }

  return comms;
}

static auto nsq_reverse_comms(std::vector<GhostCommunication> comms) {
  for (auto &c : comms) {
    c.type = GHOST_RDCE;
    if (c.comm.rank() != c.send_to)
      std::swap(c.send_lists, c.recv_lists);
    else
      c.recv_lists = c.send_lists;
  }

  return comms;
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
  cell_structure.m_local_cells.resize(1);
  cell_structure.m_local_cells[0] = local;

  cell_structure.m_ghost_cells.resize(n_nodes - 1);
  int c = 0;
  for (int n = 0; n < n_nodes; n++)
    if (n != this_node)
      cell_structure.m_ghost_cells[c++] = &cells[n];

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
  auto comms_exchange = nsq_prepare_comm(comm_cart);
  auto comms_collect = nsq_reverse_comms(comms_exchange);

  cell_structure.exchange_ghosts_comm =
      GhostCommunicator{std::move(comms_exchange)};
  cell_structure.collect_ghost_force_comm =
      GhostCommunicator{std::move(comms_collect)};

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
  displaced_parts->resize(0);

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
