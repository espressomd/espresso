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

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/range/numeric.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <vector>

/** Tag for ghosts communications. */
#define REQ_GHOST_SEND 100

/**
 * Class that stores marshalled data for ghost communications.
 * To store and retrieve data, use the adapter classes below.
 */
class CommBuf {
public:
  /** Returns a pointer to the non-bond storage.
   */
  char *data() { return buf.data(); }
  const char *data() const { return buf.data(); }

  /** Returns the number of elements in the non-bond storage.
   */
  size_t size() const { return buf.size(); }

  /** Resizes the underlying storage s.t. the object is capable
   * of holding "new_size" chars.
   * @param new_size new size
   */
  void resize(size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage.
   */
  auto &bonds() { return bondbuf; }
  const auto &bonds() const { return bondbuf; }

private:
  std::vector<char> buf;     ///< Buffer for everything but bonds
  std::vector<char> bondbuf; ///< Buffer for bond lists
};

static size_t calc_transmit_size(unsigned data_parts) {
  size_t size = {};
  if (data_parts & GHOSTTRANS_PROPRTS) {
    size += Utils::MemcpyOArchive::packing_size<ParticleProperties>();
  }
  if (data_parts & GHOSTTRANS_POSITION)
    size += Utils::MemcpyOArchive::packing_size<ParticlePosition>();
  if (data_parts & GHOSTTRANS_MOMENTUM)
    size += Utils::MemcpyOArchive::packing_size<ParticleMomentum>();
  if (data_parts & GHOSTTRANS_FORCE)
    size += Utils::MemcpyOArchive::packing_size<ParticleForce>();

  return size;
}

static size_t calc_transmit_size(GhostCommunication &ghost_comm,
                                 unsigned int data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return sizeof(int) * ghost_comm.part_lists.size();

  auto const n_part = boost::accumulate(
      ghost_comm.part_lists, 0ul,
      [](size_t sum, auto part_list) { return sum + part_list->size(); });

  return n_part * calc_transmit_size(data_parts);
}

static void prepare_send_buffer(CommBuf &send_buffer,
                                GhostCommunication &ghost_comm,
                                unsigned int data_parts) {
  /* reallocate send buffer */
  send_buffer.resize(calc_transmit_size(ghost_comm, data_parts));
  send_buffer.bonds().clear();

  auto archiver = Utils::MemcpyOArchive{Utils::make_span(send_buffer)};

  /* Construct archive that pushes back to the bond buffer */
  namespace io = boost::iostreams;
  io::stream<io::back_insert_device<std::vector<char>>> os{
      io::back_inserter(send_buffer.bonds())};
  boost::archive::binary_oarchive bond_archiver{os};

  /* put in data */
  for (auto part_list : ghost_comm.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np = part_list->size();
      archiver << np;
    } else {
      for (Particle &part : *part_list) {
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
          bond_archiver << part.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_written() == send_buffer.size());
}

static void prepare_ghost_cell(ParticleList *cell, int size) {
  /* Adapt size */
  cell->resize(size);

  /* Mark particles as ghosts */
  for (auto &p : *cell) {
    p.l.ghost = true;
  }
}

static void prepare_recv_buffer(CommBuf &recv_buffer,
                                GhostCommunication &ghost_comm,
                                unsigned int data_parts) {
  /* reallocate recv buffer */
  recv_buffer.resize(calc_transmit_size(ghost_comm, data_parts));
  /* clear bond buffer */
  recv_buffer.bonds().clear();
}

static void put_recv_buffer(CommBuf &recv_buffer,
                            GhostCommunication &ghost_comm,
                            unsigned int data_parts) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{Utils::make_span(recv_buffer)};

  if (data_parts & GHOSTTRANS_PARTNUM) {
    for (auto part_list : ghost_comm.part_lists) {
      int np;
      archiver >> np;
      prepare_ghost_cell(part_list, np);
    }
  } else {
    for (auto part_list : ghost_comm.part_lists) {
      for (Particle &part : *part_list) {
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
      }
    }
    if (data_parts & GHOSTTRANS_BONDS) {
      namespace io = boost::iostreams;
      io::stream<io::array_source> bond_stream(io::array_source{
          recv_buffer.bonds().data(), recv_buffer.bonds().size()});
      boost::archive::binary_iarchive bond_archive(bond_stream);

      for (auto part_list : ghost_comm.part_lists) {
        for (Particle &part : *part_list) {
          bond_archive >> part.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_read() == recv_buffer.size());

  recv_buffer.bonds().clear();
}

static void add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                        GhostCommunication &ghost_comm) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{Utils::make_span(recv_buffer)};
  for (auto &part_list : ghost_comm.part_lists) {
    for (Particle &part : *part_list) {
      ParticleForce pf;
      archiver >> pf;
      part.f += pf;
    }
  }
}

static void cell_cell_transfer(GhostCommunication &ghost_comm,
                               unsigned int data_parts) {
  /* transfer data */
  auto const offset = ghost_comm.part_lists.size() / 2;
  for (size_t pl = 0; pl < offset; pl++) {
    const auto *src_list = ghost_comm.part_lists[pl];
    auto *dst_list = ghost_comm.part_lists[pl + offset];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      prepare_ghost_cell(dst_list, src_list->size());
    } else {
      auto const &src_part = *src_list;
      auto &dst_part = *dst_list;
      assert(src_part.size() == dst_part.size());

      for (size_t i = 0; i < src_part.size(); i++) {
        auto const &part1 = src_part.begin()[i];
        auto &part2 = dst_part.begin()[i];

        if (data_parts & GHOSTTRANS_PROPRTS) {
          part2.p = part1.p;
        }
        if (data_parts & GHOSTTRANS_BONDS) {
          part2.bonds() = part1.bonds();
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

static bool is_send_op(int comm_type, int node, int this_node) {
  return ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
          (comm_type == GHOST_BCST && node == this_node));
}

static bool is_recv_op(int comm_type, int node, int this_node) {
  return ((comm_type == GHOST_RECV) ||
          (comm_type == GHOST_BCST && node != this_node) ||
          (comm_type == GHOST_RDCE && node == this_node));
}

static bool is_prefetchable(GhostCommunication const &ghost_comm,
                            int this_node) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const prefetch = ghost_comm.type & GHOST_PREFETCH;
  int const node = ghost_comm.node;
  return is_send_op(comm_type, node, this_node) && prefetch;
}

static bool is_poststorable(GhostCommunication const &ghost_comm,
                            int this_node) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const poststore = ghost_comm.type & GHOST_PSTSTORE;
  int const node = ghost_comm.node;
  return is_recv_op(comm_type, node, this_node) && poststore;
}

void ghost_communicator(GhostCommunicator *gcr, unsigned int data_parts) {
  if (GHOSTTRANS_NONE == data_parts)
    return;

  static CommBuf send_buffer, recv_buffer;

  auto const &comm = gcr->mpi_comm;

  for (auto it = gcr->communications.begin(); it != gcr->communications.end();
       ++it) {
    GhostCommunication &ghost_comm = *it;
    int const comm_type = ghost_comm.type & GHOST_JOBMASK;

    if (comm_type == GHOST_LOCL) {
      cell_cell_transfer(ghost_comm, data_parts);
      continue;
    }

    int const prefetch = ghost_comm.type & GHOST_PREFETCH;
    int const poststore = ghost_comm.type & GHOST_PSTSTORE;
    int const node = ghost_comm.node;

    /* prepare send buffer if necessary */
    if (is_send_op(comm_type, node, comm.rank())) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!prefetch) {
        prepare_send_buffer(send_buffer, ghost_comm, data_parts);
      }
      // Check prefetched send buffers (must also hold for buffers allocated
      // in the previous lines.)
      assert(send_buffer.size() == calc_transmit_size(ghost_comm, data_parts));
    } else if (prefetch) {
      /* we do not send this time, let's look for a prefetch */
      auto prefetch_ghost_comm = std::find_if(
          std::next(it), gcr->communications.end(),
          [this_node = comm.rank()](GhostCommunication const &ghost_comm) {
            return is_prefetchable(ghost_comm, this_node);
          });

      if (prefetch_ghost_comm != gcr->communications.end())
        prepare_send_buffer(send_buffer, *prefetch_ghost_comm, data_parts);
    }

    /* recv buffer for recv and multinode operations to this node */
    if (is_recv_op(comm_type, node, comm.rank()))
      prepare_recv_buffer(recv_buffer, ghost_comm, data_parts);

    /* transfer data */
    // Use two send/recvs in order to avoid having to serialize CommBuf
    // (which consists of already serialized data).
    switch (comm_type) {
    case GHOST_RECV:
      comm.recv(node, REQ_GHOST_SEND, recv_buffer.data(), recv_buffer.size());
      comm.recv(node, REQ_GHOST_SEND, recv_buffer.bonds());
      break;
    case GHOST_SEND:
      comm.send(node, REQ_GHOST_SEND, send_buffer.data(), send_buffer.size());
      comm.send(node, REQ_GHOST_SEND, send_buffer.bonds());
      break;
    case GHOST_BCST:
      if (node == comm.rank()) {
        boost::mpi::broadcast(comm, send_buffer.data(), send_buffer.size(),
                              node);
        boost::mpi::broadcast(comm, send_buffer.bonds(), node);
      } else {
        boost::mpi::broadcast(comm, recv_buffer.data(), recv_buffer.size(),
                              node);
        boost::mpi::broadcast(comm, recv_buffer.bonds(), node);
      }
      break;
    case GHOST_RDCE:
      if (node == comm.rank())
        boost::mpi::reduce(
            comm, reinterpret_cast<double *>(send_buffer.data()),
            static_cast<int>(send_buffer.size() / sizeof(double)),
            reinterpret_cast<double *>(recv_buffer.data()), std::plus<double>{},
            node);
      else
        boost::mpi::reduce(
            comm, reinterpret_cast<double *>(send_buffer.data()),
            static_cast<int>(send_buffer.size() / sizeof(double)),
            std::plus<double>{}, node);
      break;
    }

    // recv op; write back data directly, if no PSTSTORE delay is requested.
    if (is_recv_op(comm_type, node, comm.rank())) {
      if (!poststore) {
        /* forces have to be added, the rest overwritten. Exception is RDCE,
         * where the addition is integrated into the communication. */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, ghost_comm);
        else
          put_recv_buffer(recv_buffer, ghost_comm, data_parts);
      }
    } else if (poststore) {
      /* send op; write back delayed data from last recv, when this was a
       * prefetch send. */
      /* find previous action where we recv and which has PSTSTORE set */
      auto poststore_ghost_comm = std::find_if(
          std::make_reverse_iterator(it), gcr->communications.rend(),
          [this_node = comm.rank()](GhostCommunication const &ghost_comm) {
            return is_poststorable(ghost_comm, this_node);
          });

      if (poststore_ghost_comm != gcr->communications.rend()) {
        assert(recv_buffer.size() ==
               calc_transmit_size(*poststore_ghost_comm, data_parts));
        /* as above */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, *poststore_ghost_comm);
        else
          put_recv_buffer(recv_buffer, *poststore_ghost_comm, data_parts);
      }
    }
  }
}
