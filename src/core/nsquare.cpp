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
#include "cells.hpp"
#include "ghosts.hpp"

#include <boost/mpi/collectives/all_to_all.hpp>
#include <boost/serialization/vector.hpp>

struct AtomDecomposition {
  /** The local cell */
  Cell *local;

  boost::mpi::communicator comm;

  /**
   * @brief Determine which rank owns a particle id.
   */
  int id_to_rank(int id) const { return id % comm.size(); }

  /**
   * @brief Determine which cell a particle id belongs to.
   *
   * Since there is only one cell this is trivial.
   *
   * @param id Particle id to find cell for.
   * @return Pointer to cell or nullptr if not local.
   */
  Cell *id_to_cell(int id) {
    return (id_to_rank(id) == comm.rank()) ? local : nullptr;
  }
};

AtomDecomposition ad;

static GhostCommunicator nsq_prepare_comm() {
  /* no need for comm for only 1 node */
  if (ad.comm.size() == 1) {
    return GhostCommunicator{ad.comm, 0};
  }

  auto comm = GhostCommunicator{ad.comm, static_cast<size_t>(ad.comm.size())};
  /* every node has its dedicated comm step */
  for (int n = 0; n < ad.comm.size(); n++) {
    comm.communications[n].part_lists.resize(1);
    comm.communications[n].part_lists[0] = &(cells[n].particles());
    comm.communications[n].node = n;
  }

  return comm;
}

void nsq_topology_init(const boost::mpi::communicator &comm) {
  ad.comm = comm;

  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.particle_to_cell = [](const Particle &p) {
    return ad.id_to_cell(p.identity());
  };

  /* This cell system supports any range. */
  cell_structure.max_range =
      Utils::Vector3d::broadcast(std::numeric_limits<double>::infinity());
  cells.resize(ad.comm.size());

  /* mark cells */
  ad.local = &cells[ad.comm.rank()];
  cell_structure.m_local_cells.resize(1);
  cell_structure.m_local_cells[0] = ad.local;

  cell_structure.m_ghost_cells.resize(ad.comm.size() - 1);
  int c = 0;
  for (int n = 0; n < ad.comm.size(); n++)
    if (n != ad.comm.rank())
      cell_structure.m_ghost_cells[c++] = &cells[n];

  std::vector<Cell *> red_neighbors;
  std::vector<Cell *> black_neighbors;

  /* distribute force calculation work  */
  for (int n = 0; n < ad.comm.size(); n++) {
    auto const diff = n - ad.comm.rank();
    /* simple load balancing formula. Basically diff % n, where n >=
       ad.comm.size(), n odd. The node itself is also left out, as it is
       treated differently */
    if (diff == 0) {
      continue;
    }

    if (((diff > 0 && (diff % 2) == 0) || (diff < 0 && ((-diff) % 2) == 1))) {
      red_neighbors.push_back(&cells.at(n));
    } else {
      black_neighbors.push_back(&cells.at(n));
    }
  }

  ad.local->m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);

  /* create communicators */
  cell_structure.exchange_ghosts_comm = nsq_prepare_comm();
  cell_structure.collect_ghost_force_comm = nsq_prepare_comm();

  /* here we just decide what to transfer where */
  if (ad.comm.size() > 1) {
    for (int n = 0; n < ad.comm.size(); n++) {
      /* use the prefetched send buffers. Node 0 transmits first and never
       * prefetches. */
      if (ad.comm.rank() == 0 || ad.comm.rank() != n) {
        cell_structure.exchange_ghosts_comm.communications[n].type = GHOST_BCST;
      } else {
        cell_structure.exchange_ghosts_comm.communications[n].type =
            GHOST_BCST | GHOST_PREFETCH;
      }
      cell_structure.collect_ghost_force_comm.communications[n].type =
          GHOST_RDCE;
    }
    /* first round: all nodes except the first one prefetch their send data */
    if (ad.comm.rank() != 0) {
      cell_structure.exchange_ghosts_comm.communications[0].type |=
          GHOST_PREFETCH;
    }
  }
}

void nsq_exchange_particles(int global_flag, ParticleList *displaced_parts,
                            std::vector<Cell *> &modified_cells) {
  if (not global_flag) {
    assert(displaced_parts->empty());
    return;
  }

  /* Sort displaced particles by the node they belong to. */
  std::vector<std::vector<Particle>> send_buf(ad.comm.size());
  for (auto &p : *displaced_parts) {
    auto const target_node = ad.id_to_rank(p.identity());
    send_buf.at(target_node).emplace_back(std::move(p));
  }
  displaced_parts->clear();

  /* Exchange particles */
  std::vector<std::vector<Particle>> recv_buf(ad.comm.size());
  boost::mpi::all_to_all(ad.comm, send_buf, recv_buf);

  if (std::any_of(recv_buf.begin(), recv_buf.end(),
                  [](auto const &buf) { return not buf.empty(); })) {
    modified_cells.push_back(ad.local);
  }

  /* Add new particles belonging to this node */
  for (auto &parts : recv_buf) {
    for (auto &p : parts) {
      ad.local->particles().insert(std::move(p));
    }
  }
}
