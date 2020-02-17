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
 *  Ghost particles and particle exchange.
 *
 *  For more information on ghosts,
 *  see \ref ghosts.hpp "ghosts.hpp"
 *
 * Note on variable naming:
 * - a "GhostCommunicator" is always named "gcr",
 * - a "GhostCommunication" is always named "ghost_comm".
 */
#include "ghosts.hpp"

#include "Particle.hpp"

#include <utils/Span.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/numeric.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <vector>

/** Tag for ghosts communications. */
#define REQ_GHOST_SEND 100

static size_t calc_transmit_size(unsigned data_parts) {
  size_t size = {};
  if (data_parts & GHOSTTRANS_PROPRTS) {
    size += Utils::MemcpyOArchive::packing_size<ParticleProperties>();
  }
  if (data_parts & GHOSTTRANS_BONDS) {
    size += Utils::MemcpyOArchive::packing_size<int>();
  }
  if (data_parts & GHOSTTRANS_POSITION)
    size += Utils::MemcpyOArchive::packing_size<ParticlePosition>();
  if (data_parts & GHOSTTRANS_MOMENTUM)
    size += Utils::MemcpyOArchive::packing_size<ParticleMomentum>();
  if (data_parts & GHOSTTRANS_FORCE)
    size += Utils::MemcpyOArchive::packing_size<ParticleForce>();

  return size;
}

static size_t calc_transmit_size(std::vector<ParticleList *> const &part_lists,
                                 unsigned int data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return sizeof(int) * part_lists.size();

  auto const n_part =
      boost::accumulate(part_lists, 0ul, [](size_t sum, auto part_list) {
        return sum + part_list->n;
      });
  return n_part * calc_transmit_size(data_parts);
}

static void prepare_send_buffer(const GhostCommunication &ghost_comm,
                                unsigned int data_parts) {
  /* reallocate send buffer */
  CommBuf &send_buffer = ghost_comm.send_buffer;
  send_buffer.resize(calc_transmit_size(ghost_comm.send_lists, data_parts));
  send_buffer.bonds().clear();

  auto archiver = Utils::MemcpyOArchive{Utils::make_span(send_buffer)};
  auto bond_buffer = std::back_inserter(send_buffer.bonds());

  /* put in data */
  for (auto part_list : ghost_comm.send_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np = part_list->n;
      archiver << np;
    } else {
      for (Particle &part : part_list->particles()) {
        if (data_parts & GHOSTTRANS_PROPRTS) {
          archiver << part.p;
        }
        if (data_parts & GHOSTTRANS_POSITION) {
          /* ok, this is not nice, but perhaps fast */
          auto pp = part.r;
          pp.p += ghost_comm.shift;
          archiver << pp;
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          archiver << part.m;
        }
        if (data_parts & GHOSTTRANS_FORCE) {
          archiver << part.f;
        }
        if (data_parts & GHOSTTRANS_BONDS) {
          archiver << part.bl.n;
          boost::copy(part.bl, bond_buffer);
        }
      }
    }
  }

  assert(archiver.bytes_written() == send_buffer.size());
}

static void prepare_ghost_cell(ParticleList *cell, int size) {
  using Utils::make_span;
  auto const old_cap = cell->capacity();

  /* reset excess particles */
  if (size < cell->capacity()) {
    for (auto &p :
         make_span<Particle>(cell->part + size, cell->capacity() - size)) {
      p = Particle{};
      p.l.ghost = true;
    }
  }

  /* Adapt size */
  cell->resize(size);

  /* initialize new particles */
  if (old_cap < cell->capacity()) {
    auto new_parts =
        make_span(cell->part + old_cap, cell->capacity() - old_cap);
    std::uninitialized_fill(new_parts.begin(), new_parts.end(), Particle{});
    for (auto &p : new_parts) {
      p.l.ghost = true;
    }
  }
}

static void prepare_recv_buffer(const GhostCommunication &ghost_comm,
                                unsigned int data_parts) {
  /* reallocate recv buffer */
  ghost_comm.recv_buffer.resize(
      calc_transmit_size(ghost_comm.recv_lists, data_parts));
}

static void put_recv_buffer(const GhostCommunication &ghost_comm,
                            unsigned int data_parts) {
  /* put back data */
  CommBuf &recv_buffer = ghost_comm.recv_buffer;
  auto archiver = Utils::MemcpyIArchive{Utils::make_span(recv_buffer)};
  auto bond_buffer = recv_buffer.bonds().begin();

  for (auto part_list : ghost_comm.recv_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np;
      archiver >> np;
      prepare_ghost_cell(part_list, np);
    } else {
      for (Particle &part : part_list->particles()) {
        if (data_parts & GHOSTTRANS_PROPRTS) {
          archiver >> part.p;
        }
        if (data_parts & GHOSTTRANS_POSITION) {
          archiver >> part.r;
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          archiver >> part.m;
        }
        if (data_parts & GHOSTTRANS_FORCE) {
          archiver >> part.f;
        }
        if (data_parts & GHOSTTRANS_BONDS) {
          decltype(part.bl.n) n_bonds;
          archiver >> n_bonds;
          part.bl.resize(n_bonds);
          std::copy_n(bond_buffer, n_bonds, part.bl.begin());
          bond_buffer += n_bonds;
        }
      }
    }
  }

  assert(archiver.bytes_read() == recv_buffer.size());

  recv_buffer.bonds().clear();
}

static void add_forces_from_recv_buffer(const GhostCommunication &ghost_comm) {
  /* put back data */
  auto archiver =
      Utils::MemcpyIArchive{Utils::make_span(ghost_comm.recv_buffer)};
  for (auto &part_list : ghost_comm.recv_lists) {
    for (Particle &part : part_list->particles()) {
      ParticleForce pf;
      archiver >> pf;
      part.f += pf;
    }
  }
}

static void cell_cell_transfer(const GhostCommunication &ghost_comm,
                               unsigned int data_parts) {
  assert(ghost_comm.send_lists.size() == ghost_comm.recv_lists.size());
  /* transfer data */
  for (int pl = 0; pl < ghost_comm.send_lists.size(); pl++) {
    const ParticleList *src_list = ghost_comm.send_lists[pl];
    ParticleList *dst_list = ghost_comm.recv_lists[pl];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      prepare_ghost_cell(dst_list, src_list->n);
    } else {
      int const np = src_list->n;
      for (int p = 0; p < np; p++) {
        Particle const &part1 = src_list->part[p];
        Particle &part2 = dst_list->part[p];
        if (data_parts & GHOSTTRANS_PROPRTS) {
          part2.p = part1.p;
        }
        if (data_parts & GHOSTTRANS_BONDS) {
          part2.bl = part1.bl;
        }
        if (data_parts & GHOSTTRANS_POSITION) {
          /* ok, this is not nice, but perhaps fast */
          part2.r = part1.r;
          part2.r.p += ghost_comm.shift;
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          part2.m = part1.m;
        }
        if (data_parts & GHOSTTRANS_FORCE)
          part2.f += part1.f;
      }
    }
  }
}

void ghost_communicator(const GhostCommunicator *gcr, unsigned int data_parts) {
  for (auto const &ghost_comm : gcr->communications()) {
    int const comm_type = ghost_comm.type;

    if (comm_type == GHOST_LOCL) {
      cell_cell_transfer(ghost_comm, data_parts);
      continue;
    }

    /* prepare send buffer if necessary */
    prepare_send_buffer(ghost_comm, data_parts);

    /* recv buffer for recv and multinode operations to this node */
    prepare_recv_buffer(ghost_comm, data_parts);

    auto const send_to = ghost_comm.send_to;
    auto const recv_from = ghost_comm.recv_from;

    CommBuf &send_buffer = ghost_comm.send_buffer;
    CommBuf &recv_buffer = ghost_comm.recv_buffer;

    /* transfer data */
    // Use two send/recvs in order to avoid, having to serialize CommBuf
    // (which consists of already serialized data).
    switch (ghost_comm.type) {
    case GHOST_RECV:
      ghost_comm.comm.recv(recv_from, REQ_GHOST_SEND, recv_buffer.data(),
                           recv_buffer.size());
      ghost_comm.comm.recv(recv_from, REQ_GHOST_SEND, recv_buffer.bonds());
      break;
    case GHOST_SEND:
      ghost_comm.comm.send(send_to, REQ_GHOST_SEND, send_buffer.data(),
                           send_buffer.size());
      ghost_comm.comm.send(send_to, REQ_GHOST_SEND, send_buffer.bonds());
      break;
    case GHOST_BCST:
      if (send_to == ghost_comm.comm.rank()) {
        boost::mpi::broadcast(ghost_comm.comm, send_buffer.data(),
                              send_buffer.size(), send_to);
        boost::mpi::broadcast(ghost_comm.comm, send_buffer.bonds(), send_to);
      } else {
        boost::mpi::broadcast(ghost_comm.comm, recv_buffer.data(),
                              recv_buffer.size(), send_to);
        boost::mpi::broadcast(ghost_comm.comm, recv_buffer.bonds(), send_to);
      }
      break;
    case GHOST_RDCE:
      if (send_to == ghost_comm.comm.rank())
        boost::mpi::reduce(ghost_comm.comm,
                           reinterpret_cast<const double *>(send_buffer.data()),
                           send_buffer.size() / sizeof(double),
                           reinterpret_cast<double *>(recv_buffer.data()),
                           std::plus<double>{}, send_to);
      else
        boost::mpi::reduce(ghost_comm.comm,
                           reinterpret_cast<const double *>(send_buffer.data()),
                           send_buffer.size() / sizeof(double),
                           std::plus<double>{}, send_to);
      break;
    }

    /* forces have to be added, the rest overwritten. Exception is RDCE,
     * where the addition is integrated into the communication. */
    if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
      add_forces_from_recv_buffer(ghost_comm);
    else
      put_recv_buffer(ghost_comm, data_parts);
  }
}
