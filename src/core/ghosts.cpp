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
#include "communication.hpp"
#include "errorhandling.hpp"
#include "particle_data.hpp"

#include <boost/serialization/vector.hpp>
#include <mpi.h>
#include <utils/Span.hpp>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <type_traits>
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

  /** Returns the number of elements in the non-bond storage.
   */
  size_t size() { return buf.size(); }

  /** Resizes the underlying storage s.t. the object is capable
   * of holding "new_size" chars.
   * @param new_size new size
   */
  void resize(size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage.
   */
  std::vector<int> &bonds() { return bondbuf; }

private:
  std::vector<char> buf;    //< Buffer for everything but bonds
  std::vector<int> bondbuf; //< Buffer for bond lists
};

/**
 * Adapter to CommBuf that allows putting in and getting out data.
 * Note: The underlying buffer (span) must be large enough. This class
 * does *not* resize the buffer.
 */
struct Archiver {
private:
  Utils::Span<char> buf;
  char *insert;

public:
  explicit Archiver(Utils::Span<char> buf) : buf(buf), insert(buf.data()) {}
  ~Archiver() { assert(insert - buf.data() == buf.size()); }

  template <typename T, typename = typename std::enable_if<
                            std::is_trivially_copyable<T>::value>::type>
  void operator<<(const T &value) {
    std::memcpy(insert, &value, sizeof(T));
    insert += sizeof(T);
  }

  template <typename T, typename = typename std::enable_if<
                            std::is_trivially_copyable<T>::value>::type>
  void operator>>(T &value) {
    std::memcpy(&value, insert, sizeof(T));
    insert += sizeof(T);
  }
};

/**
 * Adapter to CommBuf that allows putting in and getting back bond data.
 * Note, this class inserts the data handed to operator<< at the end
 * of the underlying buffer. It does resize the buffer.
 */
struct BondArchiver {
private:
  /* underlying buffer type */
  using buffer_type = std::vector<int>;
  buffer_type &bondbuf;

  /* iterator into the underlying buffer */
  using const_iterator = buffer_type::const_iterator;
  const_iterator bond_retrieve;

  // Need to save these because send buffer will invalidate cb.bondbuffer.
  const const_iterator initial_begin;
  const size_t initial_size;

public:
  explicit BondArchiver(buffer_type &bondbuf)
      : bondbuf(bondbuf), bond_retrieve(bondbuf.cbegin()),
        initial_begin(bondbuf.cbegin()), initial_size(bondbuf.size()) {}

  ~BondArchiver() { assert(bond_retrieve == initial_begin + initial_size); }

  template <typename T> inline void operator<<(const Utils::Span<T> data) {
    bondbuf.insert(bondbuf.end(), data.cbegin(), data.cend());
  }

  template <typename T> inline void operator>>(Utils::Span<T> data) {
    std::copy_n(bond_retrieve, data.size(), data.begin());
    bond_retrieve += data.size();
  }
};

/** whether the ghosts should have velocity information, e.g. for DPD or RATTLE.
 *  You need this whenever you need the relative velocity of two particles.
 *  NO CHANGES OF THIS VALUE OUTSIDE OF \ref on_ghost_flags_change !!!!
 */
bool ghosts_have_v = false;
bool ghosts_have_bonds = false;

void prepare_comm(GhostCommunicator *gcr, int data_parts, int num) {
  assert(gcr);
  gcr->data_parts = data_parts;
  gcr->num = num;
  gcr->comm.resize(num);
  for (auto &ghost_comm : gcr->comm)
    ghost_comm.shift.fill(0.0);
}

void free_comm(GhostCommunicator *gcr) {
  // Invalidate the elements in all "part_lists" of all GhostCommunications.
  for (auto &ghost_comm : gcr->comm)
    ghost_comm.part_lists.clear();
}

static int calc_transmit_size(GhostCommunication &ghost_comm, int data_parts) {
  int n_buffer_new;

  if (data_parts & GHOSTTRANS_PARTNUM)
    n_buffer_new = sizeof(int) * ghost_comm.part_lists.size();
  else {
    n_buffer_new = 0;
    if (data_parts & GHOSTTRANS_PROPRTS) {
      n_buffer_new += sizeof(ParticleProperties);
      // sending size of bond lists
      if (ghosts_have_bonds) {
        n_buffer_new += sizeof(int);
      }
    }
    if (data_parts & GHOSTTRANS_POSITION)
      n_buffer_new += sizeof(ParticlePosition);
    if (data_parts & GHOSTTRANS_MOMENTUM)
      n_buffer_new += sizeof(ParticleMomentum);
    if (data_parts & GHOSTTRANS_FORCE)
      n_buffer_new += sizeof(ParticleForce);

#ifdef ENGINE
    if (data_parts & GHOSTTRANS_SWIMMING)
      n_buffer_new += sizeof(ParticleParametersSwimming);
#endif
    int count = 0;
    for (auto const &pl : ghost_comm.part_lists)
      count += pl->n;
    n_buffer_new *= count;
  }
  return n_buffer_new;
}

static void prepare_send_buffer(CommBuf &send_buffer,
                                GhostCommunication &ghost_comm,
                                int data_parts) {
  /* reallocate send buffer */
  send_buffer.resize(calc_transmit_size(ghost_comm, data_parts));
  send_buffer.bonds().clear();

  auto archiver = Archiver{Utils::make_span(send_buffer)};
  auto bond_archiver = BondArchiver{send_buffer.bonds()};

  /* put in data */
  for (auto part_list : ghost_comm.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np = part_list->n;
      archiver << np;
    } else {
      for (Particle const &part : part_list->particles()) {
        if (data_parts & GHOSTTRANS_PROPRTS) {
          archiver << part.p;
          if (ghosts_have_bonds) {
            archiver << static_cast<int>(part.bl.n);
            bond_archiver << Utils::make_const_span(part.bl);
          }
        }
        if (data_parts & GHOSTTRANS_POSSHFTD) {
          /* ok, this is not nice, but perhaps fast */
          auto pp = part.r;
          pp.p += ghost_comm.shift;
          archiver << pp;
        } else if (data_parts & GHOSTTRANS_POSITION) {
          archiver << part.r;
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          archiver << part.m;
        }
        if (data_parts & GHOSTTRANS_FORCE) {
          archiver << part.f;
        }

#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING) {
          archiver << part.swim;
        }
#endif
      }
    }
  }
}

static void prepare_ghost_cell(Cell *cell, int size) {
  using Utils::make_span;
  auto const old_cap = cell->max;

  /* reset excess particles */
  if (size < cell->max) {
    for (auto &p : make_span<Particle>(cell->part + size, cell->max - size)) {
      p = Particle{};
      p.l.ghost = true;
    }
  }

  /* Adapt size */
  cell->resize(size);

  /* initialize new particles */
  if (old_cap < cell->max) {
    auto new_parts = make_span(cell->part + old_cap, cell->max - old_cap);
    std::uninitialized_fill(new_parts.begin(), new_parts.end(), Particle{});
    for (auto &p : new_parts) {
      p.l.ghost = true;
    }
  }
}

static void prepare_recv_buffer(CommBuf &recv_buffer,
                                GhostCommunication &ghost_comm,
                                int data_parts) {
  /* reallocate recv buffer */
  recv_buffer.resize(calc_transmit_size(ghost_comm, data_parts));
}

static void put_recv_buffer(CommBuf &recv_buffer,
                            GhostCommunication &ghost_comm, int data_parts) {
  /* put back data */
  auto archiver = Archiver{Utils::make_span(recv_buffer)};
  auto bond_archiver = BondArchiver{recv_buffer.bonds()};

  for (auto part_list : ghost_comm.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np;
      archiver >> np;
      prepare_ghost_cell(part_list, np);
    } else {
      for (Particle &part : part_list->particles()) {
        if (data_parts & GHOSTTRANS_PROPRTS) {
          archiver >> part.p;
          if (ghosts_have_bonds) {
            int n_bonds;
            archiver >> n_bonds;
            part.bl.resize(n_bonds);
            bond_archiver >> Utils::make_span(part.bl);
          }
          if (local_particles[part.p.identity] == nullptr) {
            local_particles[part.p.identity] = &part;
          }
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

#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING) {
          archiver >> part.swim;
        }
#endif
      }
    }
  }

  recv_buffer.bonds().clear();
}

static void add_forces_from_recv_buffer(CommBuf &recv_buffer,
                                        GhostCommunication &ghost_comm) {
  /* put back data */
  auto archiver = Archiver{Utils::make_span(recv_buffer)};
  for (auto &part_list : ghost_comm.part_lists) {
    for (Particle &part : part_list->particles()) {
      ParticleForce pf;
      archiver >> pf;
      part.f += pf;
    }
  }
}

static void cell_cell_transfer(GhostCommunication *ghost_comm, int data_parts) {
  /* transfer data */
  int const offset = ghost_comm->part_lists.size() / 2;
  for (int pl = 0; pl < offset; pl++) {
    Cell *src_list = ghost_comm->part_lists[pl];
    Cell *dst_list = ghost_comm->part_lists[pl + offset];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      prepare_ghost_cell(dst_list, src_list->n);
    } else {
      int const np = src_list->n;
      for (int p = 0; p < np; p++) {
        Particle const &part1 = src_list->part[p];
        Particle &part2 = dst_list->part[p];
        if (data_parts & GHOSTTRANS_PROPRTS) {
          part2.p = part1.p;
          if (ghosts_have_bonds) {
            part2.bl = part1.bl;
          }
        }
        if (data_parts & GHOSTTRANS_POSSHFTD) {
          /* ok, this is not nice, but perhaps fast */
          part2.r = part1.r;
          part2.r.p += ghost_comm->shift;
        } else if (data_parts & GHOSTTRANS_POSITION)
          part2.r = part1.r;
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          part2.m = part1.m;
        }
        if (data_parts & GHOSTTRANS_FORCE)
          part2.f += part1.f;

#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING)
          part2.swim = part1.swim;
#endif
      }
    }
  }
}

static bool is_send_op(int comm_type, int node) {
  return ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
          (comm_type == GHOST_BCST && node == this_node));
}

static bool is_recv_op(int comm_type, int node) {
  return ((comm_type == GHOST_RECV) ||
          (comm_type == GHOST_BCST && node != this_node) ||
          (comm_type == GHOST_RDCE && node == this_node));
}

static bool is_prefetchable(GhostCommunication const &ghost_comm) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const prefetch = ghost_comm.type & GHOST_PREFETCH;
  int const node = ghost_comm.node;
  return is_send_op(comm_type, node) && prefetch;
}

static bool is_poststorable(GhostCommunication const &ghost_comm) {
  int const comm_type = ghost_comm.type & GHOST_JOBMASK;
  int const poststore = ghost_comm.type & GHOST_PSTSTORE;
  int const node = ghost_comm.node;
  return is_recv_op(comm_type, node) && poststore;
}

void ghost_communicator(GhostCommunicator *gcr) {
  ghost_communicator(gcr, gcr->data_parts);
}

void ghost_communicator(GhostCommunicator *gcr, int data_parts) {
  static CommBuf send_buffer, recv_buffer;

  /* if ghosts should have uptodate velocities, they have to be updated like
   * positions (except for shifting...) */
  if (ghosts_have_v && (data_parts & GHOSTTRANS_POSITION))
    data_parts |= GHOSTTRANS_MOMENTUM;

  for (int n = 0; n < gcr->num; n++) {
    GhostCommunication *ghost_comm = &gcr->comm[n];
    int const comm_type = ghost_comm->type & GHOST_JOBMASK;

    if (comm_type == GHOST_LOCL) {
      cell_cell_transfer(ghost_comm, data_parts);
      continue;
    }

    int const prefetch = ghost_comm->type & GHOST_PREFETCH;
    int const poststore = ghost_comm->type & GHOST_PSTSTORE;
    int const node = ghost_comm->node;

    /* prepare send buffer if necessary */
    if (is_send_op(comm_type, node)) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!prefetch)
        prepare_send_buffer(send_buffer, *ghost_comm, data_parts);
#ifdef ADDITIONAL_CHECKS
      // Check prefetched send buffers (must also hold for buffers allocated
      // in the previous lines.)
      if (send_buffer.size() != calc_transmit_size(*ghost_comm, data_parts)) {
        fprintf(stderr,
                "%d: ghost_comm transmission size and current size of "
                "cells to transmit do not match\n",
                this_node);
        errexit();
      }
#endif
    } else if (prefetch) {
      /* we do not send this time, let's look for a prefetch */
      auto prefetch_ghost_comm = std::find_if(
          std::next(gcr->comm.begin(), n + 1), gcr->comm.end(), is_prefetchable);
      if (prefetch_ghost_comm != gcr->comm.end())
        prepare_send_buffer(send_buffer, *prefetch_ghost_comm, data_parts);
    }

    /* recv buffer for recv and multinode operations to this node */
    if (is_recv_op(comm_type, node))
      prepare_recv_buffer(recv_buffer, *ghost_comm, data_parts);

    /* transfer data */
    // Use two send/recvs in order to avoid, having to serialize CommBuf
    // (which consists of already serialized data).
    switch (comm_type) {
    case GHOST_RECV:
      comm_cart.recv(node, REQ_GHOST_SEND, recv_buffer.data(),
                     recv_buffer.size());
      comm_cart.recv(node, REQ_GHOST_SEND, recv_buffer.bonds());
      break;
    case GHOST_SEND:
      comm_cart.send(node, REQ_GHOST_SEND, send_buffer.data(),
                     send_buffer.size());
      comm_cart.send(node, REQ_GHOST_SEND, send_buffer.bonds());
      break;
    case GHOST_BCST:
      if (node == this_node) {
        boost::mpi::broadcast(comm_cart, send_buffer.data(), send_buffer.size(),
                              node);
        boost::mpi::broadcast(comm_cart, send_buffer.bonds(), node);
      } else {
        boost::mpi::broadcast(comm_cart, recv_buffer.data(), recv_buffer.size(),
                              node);
        boost::mpi::broadcast(comm_cart, recv_buffer.bonds(), node);
      }
      break;
    case GHOST_RDCE:
      if (node == this_node)
        boost::mpi::reduce(comm_cart,
                           reinterpret_cast<double *>(send_buffer.data()),
                           send_buffer.size() / sizeof(double),
                           reinterpret_cast<double *>(recv_buffer.data()),
                           std::plus<double>{}, node);
      else
        boost::mpi::reduce(
            comm_cart, reinterpret_cast<double *>(send_buffer.data()),
            send_buffer.size() / sizeof(double), std::plus<double>{}, node);
      break;
    }

    // recv op; write back data directly, if no PSTSTORE delay is requested.
    if (is_recv_op(comm_type, node)) {
      if (!poststore) {
        /* forces have to be added, the rest overwritten. Exception is RDCE,
         * where the addition is integrated into the communication. */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, *ghost_comm);
        else
          put_recv_buffer(recv_buffer, *ghost_comm, data_parts);
      }
    } else if (poststore) {
      /* send op; write back delayed data from last recv, when this was a
       * prefetch send. */
      /* find previous action where we recv and which has PSTSTORE set */
      auto poststore_ghost_comm = std::find_if(
          std::make_reverse_iterator(std::next(gcr->comm.begin(), n)),
          gcr->comm.rend(), is_poststorable);

      if (poststore_ghost_comm != gcr->comm.rend()) {
#ifdef ADDITIONAL_CHECKS
        if (recv_buffer.size() !=
            calc_transmit_size(*poststore_ghost_comm, data_parts)) {
          fprintf(stderr,
                  "%d: ghost_comm transmission size and current size of "
                  "cells to transmit do not match\n",
                  this_node);
          errexit();
        }
#endif
        /* as above */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(recv_buffer, *poststore_ghost_comm);
        else
          put_recv_buffer(recv_buffer, *poststore_ghost_comm, data_parts);
      }
    }
  }
}
