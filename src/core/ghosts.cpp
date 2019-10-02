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

/** Tag for communication in ghost_comm. */
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
  std::vector<char> buf; //< Buffer for everything but bonds
  std::vector<int> bondbuf; //< Buffer for bond lists
};

/**
 * Adapter to CommBuf that allows putting in and getting out data.
 * Note: The underlying buffer (span) must be large enough. This class
 * does *not* resize the buffer.
 */
struct Archiver {
private:
  Utils::Span<char> cb;
  char *insert;

public:
  explicit Archiver(Utils::Span<char> cb) : cb(cb), insert(cb.data()) {}
  ~Archiver() { assert(insert - cb.data() == cb.size()); }

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

void prepare_comm(GhostCommunicator *comm, int data_parts, int num) {
  assert(comm);
  comm->data_parts = data_parts;
  comm->num = num;
  comm->comm.resize(num);
  for (auto &gc : comm->comm)
    gc.shift.fill(0.0);
}

void free_comm(GhostCommunicator *comm) {
  // Invalidate the elements in all "part_lists" of all GhostCommunications.
  for (auto &gc : comm->comm)
    gc.part_lists.clear();
}

static int calc_transmit_size(GhostCommunication &gc, int data_parts) {
  int n_buffer_new;

  if (data_parts & GHOSTTRANS_PARTNUM)
    n_buffer_new = sizeof(int) * gc.part_lists.size();
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
    for (auto const &pl : gc.part_lists)
      count += pl->n;
    n_buffer_new *= count;
  }
  return n_buffer_new;
}

static void prepare_send_buffer(CommBuf &s_buffer, GhostCommunication &gc,
                                int data_parts) {
  /* reallocate send buffer */
  s_buffer.resize(calc_transmit_size(gc, data_parts));
  s_buffer.bonds().clear();

  auto ar = Archiver{Utils::make_span(s_buffer)};
  auto bar = BondArchiver{s_buffer.bonds()};

  /* put in data */
  for (auto cur_list : gc.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np = cur_list->n;
      ar << np;
    } else {
      for (Particle const &pt : cur_list->particles()) {
        if (data_parts & GHOSTTRANS_PROPRTS) {
          ar << pt.p;
          if (ghosts_have_bonds) {
            ar << static_cast<int>(pt.bl.n);
            bar << Utils::make_const_span(pt.bl);
          }
        }
        if (data_parts & GHOSTTRANS_POSSHFTD) {
          /* ok, this is not nice, but perhaps fast */
          auto pp = pt.r;
          pp.p += gc.shift;
          ar << pp;
        } else if (data_parts & GHOSTTRANS_POSITION) {
          ar << pt.r;
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          ar << pt.m;
        }
        if (data_parts & GHOSTTRANS_FORCE) {
          ar << pt.f;
        }

#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING) {
          ar << pt.swim;
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

static void prepare_recv_buffer(CommBuf &r_buffer, GhostCommunication &gc,
                                int data_parts) {
  /* reallocate recv buffer */
  r_buffer.resize(calc_transmit_size(gc, data_parts));
}

static void put_recv_buffer(CommBuf &r_buffer, GhostCommunication &gc,
                            int data_parts) {
  /* put back data */
  auto ar = Archiver{Utils::make_span(r_buffer)};
  auto bar = BondArchiver{r_buffer.bonds()};

  for (auto cur_list : gc.part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      int np;
      ar >> np;
      prepare_ghost_cell(cur_list, np);
    } else {
      for (Particle &pt : cur_list->particles()) {
        if (data_parts & GHOSTTRANS_PROPRTS) {
          ar >> pt.p;
          if (ghosts_have_bonds) {
            int n_bonds;
            ar >> n_bonds;
            pt.bl.resize(n_bonds);
            bar >> Utils::make_span(pt.bl);
          }
          if (local_particles[pt.p.identity] == nullptr) {
            local_particles[pt.p.identity] = &pt;
          }
        }
        if (data_parts & GHOSTTRANS_POSITION) {
          ar >> pt.r;
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          ar >> pt.m;
        }
        if (data_parts & GHOSTTRANS_FORCE) {
          ar >> pt.f;
        }

#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING) {
          ar >> pt.swim;
        }
#endif
      }
    }
  }

  r_buffer.bonds().clear();
}

static void add_forces_from_recv_buffer(CommBuf &r_buffer,
                                        GhostCommunication &gc) {
  /* put back data */
  auto ar = Archiver{Utils::make_span(r_buffer)};
  for (auto &pl : gc.part_lists) {
    for (Particle &pt : pl->particles()) {
      ParticleForce pf;
      ar >> pf;
      pt.f += pf;
    }
  }
}

static void cell_cell_transfer(GhostCommunication *gc, int data_parts) {
  /* transfer data */
  int const offset = gc->part_lists.size() / 2;
  for (int pl = 0; pl < offset; pl++) {
    Cell *src_list = gc->part_lists[pl];
    Cell *dst_list = gc->part_lists[pl + offset];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      prepare_ghost_cell(dst_list, src_list->n);
    } else {
      int const np = src_list->n;
      for (int p = 0; p < np; p++) {
        Particle const &pt1 = src_list->part[p];
        Particle &pt2 = dst_list->part[p];
        if (data_parts & GHOSTTRANS_PROPRTS) {
          pt2.p = pt1.p;
          if (ghosts_have_bonds) {
            pt2.bl = pt1.bl;
          }
        }
        if (data_parts & GHOSTTRANS_POSSHFTD) {
          /* ok, this is not nice, but perhaps fast */
          pt2.r = pt1.r;
          pt2.r.p += gc->shift;
        } else if (data_parts & GHOSTTRANS_POSITION)
          pt2.r = pt1.r;
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          pt2.m = pt1.m;
        }
        if (data_parts & GHOSTTRANS_FORCE)
          pt2.f += pt1.f;

#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING)
          pt2.swim = pt1.swim;
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

static bool is_prefetchable(GhostCommunication const &gcn) {
  int const comm_type = gcn.type & GHOST_JOBMASK;
  int const prefetch = gcn.type & GHOST_PREFETCH;
  int const node = gcn.node;
  return is_send_op(comm_type, node) && prefetch;
}

static bool is_poststorable(GhostCommunication const &gcn) {
  int const comm_type = gcn.type & GHOST_JOBMASK;
  int const poststore = gcn.type & GHOST_PSTSTORE;
  int const node = gcn.node;
  return is_recv_op(comm_type, node) && poststore;
}

void ghost_communicator(GhostCommunicator *gc) {
  ghost_communicator(gc, gc->data_parts);
}

void ghost_communicator(GhostCommunicator *gc, int data_parts) {
  static CommBuf s_buffer, r_buffer;

  /* if ghosts should have uptodate velocities, they have to be updated like
   * positions (except for shifting...) */
  if (ghosts_have_v && (data_parts & GHOSTTRANS_POSITION))
    data_parts |= GHOSTTRANS_MOMENTUM;

  for (int n = 0; n < gc->num; n++) {
    GhostCommunication *gcn = &gc->comm[n];
    int const comm_type = gcn->type & GHOST_JOBMASK;

    if (comm_type == GHOST_LOCL) {
      cell_cell_transfer(gcn, data_parts);
      continue;
    }

    int const prefetch = gcn->type & GHOST_PREFETCH;
    int const poststore = gcn->type & GHOST_PSTSTORE;
    int const node = gcn->node;

    /* prepare send buffer if necessary */
    if (is_send_op(comm_type, node)) {
      /* ok, we send this step, prepare send buffer if not yet done */
      if (!prefetch)
        prepare_send_buffer(s_buffer, *gcn, data_parts);
#ifdef ADDITIONAL_CHECKS
      // Check prefetched send buffers (must also hold for buffers allocated
      // in the previous lines.)
      if (s_buffer.size() != calc_transmit_size(*gcn, data_parts)) {
        fprintf(stderr,
                "%d: ghost_comm transmission size and current size of "
                "cells to transmit do not match\n",
                this_node);
        errexit();
      }
#endif
    } else if (prefetch) {
      /* we do not send this time, let's look for a prefetch */
      auto pref_gcn = std::find_if(std::next(gc->comm.begin(), n + 1),
                                   gc->comm.end(), is_prefetchable);
      if (pref_gcn != gc->comm.end())
        prepare_send_buffer(s_buffer, *pref_gcn, data_parts);
    }

    /* recv buffer for recv and multinode operations to this node */
    if (is_recv_op(comm_type, node))
      prepare_recv_buffer(r_buffer, *gcn, data_parts);

    /* transfer data */
    // Use two send/recvs in order to avoid, having to serialize CommBuf
    // (which consists of already serialized data).
    switch (comm_type) {
    case GHOST_RECV:
      comm_cart.recv(node, REQ_GHOST_SEND, r_buffer.data(), r_buffer.size());
      comm_cart.recv(node, REQ_GHOST_SEND, r_buffer.bonds());
      break;
    case GHOST_SEND:
      comm_cart.send(node, REQ_GHOST_SEND, s_buffer.data(), s_buffer.size());
      comm_cart.send(node, REQ_GHOST_SEND, s_buffer.bonds());
      break;
    case GHOST_BCST:
      if (node == this_node) {
        boost::mpi::broadcast(comm_cart, s_buffer.data(), s_buffer.size(),
                              node);
        boost::mpi::broadcast(comm_cart, s_buffer.bonds(), node);
      } else {
        boost::mpi::broadcast(comm_cart, r_buffer.data(), r_buffer.size(),
                              node);
        boost::mpi::broadcast(comm_cart, r_buffer.bonds(), node);
      }
      break;
    case GHOST_RDCE:
      if (node == this_node)
        boost::mpi::reduce(comm_cart,
                           reinterpret_cast<double *>(s_buffer.data()),
                           s_buffer.size() / sizeof(double),
                           reinterpret_cast<double *>(r_buffer.data()),
                           std::plus<double>{}, node);
      else
        boost::mpi::reduce(
            comm_cart, reinterpret_cast<double *>(s_buffer.data()),
            s_buffer.size() / sizeof(double), std::plus<double>{}, node);
      break;
    }

    // recv op; write back data directly, if no PSTSTORE delay is requested.
    if (is_recv_op(comm_type, node)) {
      if (!poststore) {
        /* forces have to be added, the rest overwritten. Exception is RDCE,
         * where the addition is integrated into the communication. */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(r_buffer, *gcn);
        else
          put_recv_buffer(r_buffer, *gcn, data_parts);
      }
    } else if (poststore) {
      /* send op; write back delayed data from last recv, when this was a
       * prefetch send. */
      /* find previous action where we recv and which has PSTSTORE set */
      auto postst_gcn = std::find_if(
          std::make_reverse_iterator(std::next(gc->comm.begin(), n)),
          gc->comm.rend(), is_poststorable);

      if (postst_gcn != gc->comm.rend()) {
#ifdef ADDITIONAL_CHECKS
        if (r_buffer.size() != calc_transmit_size(*postst_gcn, data_parts)) {
          fprintf(stderr,
                  "%d: ghost_comm transmission size and current size of "
                  "cells to transmit do not match\n",
                  this_node);
          errexit();
        }
#endif
        /* as above */
        if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
          add_forces_from_recv_buffer(r_buffer, *postst_gcn);
        else
          put_recv_buffer(r_buffer, *postst_gcn, data_parts);
      }
    }
  }
}
