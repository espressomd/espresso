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
 *  Ghost particles and particle exchange.
 *
 *  For more information on ghosts,
 *  see \ref ghosts.hpp "ghosts.hpp"
 */
#include "ghosts.hpp"
#include "cells.hpp"
#include "debug.hpp"
#include "particle_data.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/range/numeric.hpp>
#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

/** Tag for communication in ghost_comm. */
#define REQ_GHOST_SEND 100

static size_t n_s_buffer = 0;
static size_t max_s_buffer = 0;
/** send buffer. Just grows, which should be ok */
static char *s_buffer = nullptr;

static std::vector<int> s_bondbuffer;

static size_t n_r_buffer = 0;
static size_t max_r_buffer = 0;
/** recv buffer. Just grows, which should be ok */
static char *r_buffer = nullptr;

static std::vector<int> r_bondbuffer;

/** whether the ghosts should also have velocity information, e. g. for DPD or
   RATTLE. You need this whenever you need the relative velocity of two
   particles. NO CHANGES OF THIS VALUE OUTSIDE OF \ref on_ghost_flags_change
   !!!!
*/
int ghosts_have_v = 0;
int ghosts_have_bonds = 0;

static size_t calc_size_per_part(int data_parts) {
  size_t size = 0;
  if (data_parts & GHOSTTRANS_PROPRTS) {
    size += sizeof(ParticleProperties);
    // sending size of bond lists
    if (ghosts_have_bonds) {
      size += sizeof(int);
    }
  }
  if (data_parts & GHOSTTRANS_POSITION)
    size += sizeof(ParticlePosition);
  if (data_parts & GHOSTTRANS_MOMENTUM)
    size += sizeof(ParticleMomentum);
  if (data_parts & GHOSTTRANS_FORCE)
    size += sizeof(ParticleForce);
#ifdef LB
  if (data_parts & GHOSTTRANS_COUPLING)
    size += sizeof(ParticleLatticeCoupling);
#endif
#ifdef ENGINE
  if (data_parts & GHOSTTRANS_SWIMMING)
    size += sizeof(ParticleParametersSwimming);
#endif
  return size;
}

static size_t calc_transmit_size(const std::vector<Cell *> &part_lists, int data_parts) {
  size_t n_buffer_new = 0;

  if (data_parts & GHOSTTRANS_PARTNUM) {
    n_buffer_new = sizeof(int) * part_lists.size();
  } else {
    auto const count = boost::accumulate(
            part_lists, 0, [](int sum, const Cell *c) { return sum + c->n; });

    auto const size_per_part = calc_size_per_part(data_parts);
    n_buffer_new = count * size_per_part;
  }
  // also sending length of bond buffer
  if (data_parts & GHOSTTRANS_PROPRTS)
    n_buffer_new += sizeof(int);
  return n_buffer_new;
}

static char * pack_particle(const Particle & pt, char * insert, int data_parts, const boost::optional<Vector3d> & shift) {
  if (data_parts & GHOSTTRANS_PROPRTS) {
    memcpy(insert, &pt.p, sizeof(ParticleProperties));
    insert += sizeof(ParticleProperties);
    if (ghosts_have_bonds) {
      *(int *)insert = pt.bl.n;
      insert += sizeof(int);
    }
  }
  if (data_parts & GHOSTTRANS_POSITION) {
    auto pp = new (insert) ParticlePosition(pt.r);

    if (shift) {
      pp->p += shift.get();
    }
    insert += sizeof(ParticlePosition);
  }
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    memcpy(insert, &pt.m, sizeof(ParticleMomentum));
    insert += sizeof(ParticleMomentum);
  }
  if (data_parts & GHOSTTRANS_FORCE) {
    memcpy(insert, &pt.f, sizeof(ParticleForce));
    insert += sizeof(ParticleForce);
  }
#ifdef LB
  if (data_parts & GHOSTTRANS_COUPLING) {
    memcpy(insert, &pt.lc, sizeof(ParticleLatticeCoupling));
    insert += sizeof(ParticleLatticeCoupling);
  }
#endif
#ifdef ENGINE
  if (data_parts & GHOSTTRANS_SWIMMING) {
    memcpy(insert, &pt.swim, sizeof(ParticleParametersSwimming));
    insert += sizeof(ParticleParametersSwimming);
  }
#endif
  return insert;
}

static void prepare_send_buffer(GhostCommunication *gc, int data_parts, boost::optional<Vector3d> const& shift) {
  GHOST_TRACE(fprintf(stderr, "%d: prepare sending to/bcast from %d\n",
                      this_node, gc->node));

  /* reallocate send buffer */
  n_s_buffer = calc_transmit_size(gc->part_lists, data_parts);
  if (n_s_buffer > max_s_buffer) {
    max_s_buffer = n_s_buffer;
    s_buffer = Utils::realloc(s_buffer, max_s_buffer);
  }
  GHOST_TRACE(fprintf(stderr, "%d: will send %d\n", this_node, n_s_buffer));

  s_bondbuffer.clear();

  /* put in data */
  char *insert = s_buffer;
  for (auto const &pl : gc->part_lists) {
    int np = pl->n;
    if (data_parts & GHOSTTRANS_PARTNUM) {
      *(int *)insert = np;
      insert += sizeof(int);
      GHOST_TRACE(
          fprintf(stderr, "%d: %d particles assigned\n", this_node, np));
    } else {
      for (int p = 0; p < np; p++) {
        const Particle &pt = pl->part[p];
        insert = pack_particle(pt, insert, data_parts, shift);

        if (ghosts_have_bonds and (data_parts & GHOSTTRANS_PROPRTS)) {
          if (pt.bl.n) {
            s_bondbuffer.insert(s_bondbuffer.end(), pt.bl.e,
                                pt.bl.e + pt.bl.n);
          }
        }
      }
    }
  }

  if (data_parts & GHOSTTRANS_PROPRTS) {
    GHOST_TRACE(fprintf(stderr, "%d: bond buffer size is %ld\n", this_node,
                        s_bondbuffer.size()));
    *(int *)insert = int(s_bondbuffer.size());
    insert += sizeof(int);
  }

  /* Assert that the calculated and actual buffer size match. */
  assert((insert - s_buffer) == n_s_buffer);
}

static void prepare_ghost_cell(Cell *cell, int size) {
  if (size > cell->max) {
    auto const old_cap = cell->max;
    cell->max = size;
    cell->part = Utils::realloc(cell->part, size * sizeof(Particle));

    auto p = Particle();
    p.p.identity = -1;
    p.l.ghost = 1;

    std::uninitialized_fill(cell->part + old_cap, cell->part + cell->max, p);
  }

  cell->n = size;
}

static void prepare_recv_buffer(GhostCommunication *gc, int data_parts) {
  GHOST_TRACE(
      fprintf(stderr, "%d: prepare receiving from %d\n", this_node, gc->node));
  /* reallocate recv buffer */
  n_r_buffer = calc_transmit_size(gc->part_lists, data_parts);
  if (n_r_buffer > max_r_buffer) {
    max_r_buffer = n_r_buffer;
    r_buffer = Utils::realloc(r_buffer, max_r_buffer);
  }
  GHOST_TRACE(fprintf(stderr, "%d: will get %d\n", this_node, n_r_buffer));
}

static void put_recv_buffer(GhostCommunication *gc, int data_parts) {
  /* put back data */
  char *retrieve = r_buffer;

  std::vector<int>::const_iterator bond_retrieve = r_bondbuffer.begin();

  for (auto &pl : gc->part_lists) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      GHOST_TRACE(fprintf(
          stderr, "%d: reallocating cell %p to size %d, assigned to node %d\n",
          this_node, (void *)cur_list, *(int *)retrieve, gc->node));
      prepare_ghost_cell(pl, *(int *)retrieve);
      retrieve += sizeof(int);
    } else {
      int np = pl->n;
      Particle *part = pl->part;
      for (int p = 0; p < np; p++) {
        Particle *pt = &part[p];
        if (data_parts & GHOSTTRANS_PROPRTS) {
          memcpy(&pt->p, retrieve, sizeof(ParticleProperties));
          retrieve += sizeof(ParticleProperties);
          if (ghosts_have_bonds) {
            int n_bonds;
            memcpy(&n_bonds, retrieve, sizeof(int));
            retrieve += sizeof(int);
            if (n_bonds) {
              pt->bl.resize(n_bonds);
              std::copy_n(bond_retrieve, n_bonds, pt->bl.begin());
              bond_retrieve += n_bonds;
            }
          }
          if (local_particles[pt->p.identity] == nullptr) {
            local_particles[pt->p.identity] = pt;
          }
        }
        if (data_parts & GHOSTTRANS_POSITION) {
          memcpy(&pt->r, retrieve, sizeof(ParticlePosition));
          retrieve += sizeof(ParticlePosition);
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          memcpy(&pt->m, retrieve, sizeof(ParticleMomentum));
          retrieve += sizeof(ParticleMomentum);
        }
        if (data_parts & GHOSTTRANS_FORCE) {
          memcpy(&pt->f, retrieve, sizeof(ParticleForce));
          retrieve += sizeof(ParticleForce);
        }
#ifdef LB
        if (data_parts & GHOSTTRANS_COUPLING) {
          memcpy(&pt->lc, retrieve, sizeof(ParticleLatticeCoupling));
          retrieve += sizeof(ParticleLatticeCoupling);
        }
#endif
#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING) {
          memcpy(&pt->swim, retrieve, sizeof(ParticleParametersSwimming));
          retrieve += sizeof(ParticleParametersSwimming);
        }
#endif
      }
    }
  }

  if (data_parts & GHOSTTRANS_PROPRTS) {
    // skip the final information on bonds to be sent in a second round
    retrieve += sizeof(int);
  }

  /* Assert that the actual and calculated buffer size match */
  assert((retrieve - r_buffer) == n_r_buffer);
  /* Assert that all bonds have been consumed. */
  assert(bond_retrieve == r_bondbuffer.end());
  r_bondbuffer.clear();
}

static void add_forces_from_recv_buffer(std::vector<Cell *> const& part_lists,
        Utils::Span<const ParticleForce> buffer) {
  auto it = buffer.begin();

  for (auto &pl : part_lists) {
    int np = pl->n;
    auto part = pl->part;
    for (int p = 0; p < np; p++) {
      auto &pt = part[p];
      pt.f += *it++;
    }
  }

  assert(it == buffer.end());
}

static void cell_cell_transfer(GhostCommunication *gc, int data_parts) {
  int pl, p;
  Particle *part1, *part2, *pt1, *pt2;

  GHOST_TRACE(fprintf(stderr, "%d: local_transfer: type %d data_parts %d\n",
                      this_node, gc->type, data_parts));

  /* transfer data */
  auto const offset = gc->part_lists.size() / 2;
  for (pl = 0; pl < offset; pl++) {
    Cell *src_list = gc->part_lists[pl];
    Cell *dst_list = gc->part_lists[pl + offset];

    if (data_parts & GHOSTTRANS_PARTNUM) {
      prepare_ghost_cell(dst_list, src_list->n);
    } else {
      int np = src_list->n;
      part1 = src_list->part;
      part2 = dst_list->part;
      for (p = 0; p < np; p++) {
        pt1 = &part1[p];
        pt2 = &part2[p];
        if (data_parts & GHOSTTRANS_PROPRTS) {
          pt2->p = pt1->p;
          if (ghosts_have_bonds) {
            pt2->bl = pt1->bl;
          }
        }
        if (data_parts & GHOSTTRANS_POSITION) {
          pt2->r = pt1->r;
          if (gc->shift) {
            pt2->r.p += gc->shift.get();
          }
        }
        if (data_parts & GHOSTTRANS_MOMENTUM) {
          pt2->m = pt1->m;
        }
        if (data_parts & GHOSTTRANS_FORCE)
          pt2->f += pt1->f;

#ifdef LB
        if (data_parts & GHOSTTRANS_COUPLING)
          pt2->lc = pt1->lc;
#endif
#ifdef ENGINE
        if (data_parts & GHOSTTRANS_SWIMMING)
          pt2->swim = pt1->swim;
#endif
      }
    }
  }
}

static int is_send_op(int comm_type, int node) {
  return ((comm_type == GHOST_SEND) || (comm_type == GHOST_RDCE) ||
          (comm_type == GHOST_BCST && node == this_node));
}

static int is_recv_op(int comm_type, int node) {
  return ((comm_type == GHOST_RECV) ||
          (comm_type == GHOST_BCST && node != this_node) ||
          (comm_type == GHOST_RDCE && node == this_node));
}

namespace {
template <class Iter>
constexpr std::reverse_iterator<Iter> make_reverse_iterator(Iter i) {
  return std::reverse_iterator<Iter>(i);
}
} // namespace

void ghost_communicator(GhostCommunicator &gc, int data_parts) {
  /* if ghosts should have up-to-date velocities, they have to be updated like
     positions (except for shifting...) */
  if (ghosts_have_v && (data_parts & GHOSTTRANS_POSITION))
    data_parts |= GHOSTTRANS_MOMENTUM;

  GHOST_TRACE(fprintf(stderr, "%d: ghost_comm %p, data_parts %d\n", this_node,
                      (void *)gc, data_parts));

  for (auto it = gc.comm.begin(); it != gc.comm.end(); ++it) {
    GhostCommunication *gcn = &(*it);
    int comm_type = gcn->type;
    int prefetch = gcn->prefetch;
    int poststore = gcn->poststore;
    int node = gcn->node;

    if (comm_type == GHOST_LOCL)
      cell_cell_transfer(gcn, data_parts);
    else {
      /* prepare send buffer if necessary */
      if (is_send_op(comm_type, node)) {
        /* ok, we send this step, prepare send buffer if not yet done */
        if (!prefetch)
          prepare_send_buffer(gcn, data_parts, gcn->shift);
        else {
          GHOST_TRACE(fprintf(stderr,
                              "%d: ghost_comm using prefetched data for "
                              "operation %d, sending to %d\n",
                              this_node, n, node));
          assert(n_s_buffer == calc_transmit_size(gcn->part_lists, data_parts));
        }
      } else {
        /* we do not send this time, let's look for a prefetch */
        if (prefetch) {
          /* find next action where we send and which has PREFETCH set */
          for (auto jt = std::next(it); jt != gc.comm.end(); ++jt) {
            GhostCommunication *gcn2 = &(*jt);
            int comm_type2 = gcn2->type;
            int prefetch2 = gcn2->prefetch;
            int node2 = gcn2->node;
            if (is_send_op(comm_type2, node2) && prefetch2) {
              prepare_send_buffer(gcn2, data_parts, gcn2->shift);
              break;
            }
          }
        }
      }

      /* recv buffer for recv and multinode operations to this node */
      if (is_recv_op(comm_type, node))
        prepare_recv_buffer(gcn, data_parts);

      /* transfer data */
      switch (comm_type) {
      case GHOST_RECV: {
        GHOST_TRACE(fprintf(stderr,
                            "%d: ghost_comm receive from %d (%d bytes)\n",
                            this_node, node, n_r_buffer));
        if (n_r_buffer > 0) {
          gc.mpi_comm.recv(node, REQ_GHOST_SEND, r_buffer, n_r_buffer);
        }

        if (data_parts & GHOSTTRANS_PROPRTS) {
          int n_bonds = *(int *)(r_buffer + n_r_buffer - sizeof(int));
          GHOST_TRACE(fprintf(stderr,
                              "%d: ghost_comm receive from %d (%d bonds)\n",
                              this_node, node, n_bonds));
          if (n_bonds) {
            r_bondbuffer.resize(n_bonds);
            gc.mpi_comm.recv(node, REQ_GHOST_SEND, r_bondbuffer.data(),
                             n_bonds);
          }
        }
        break;
      }
      case GHOST_SEND: {
        GHOST_TRACE(fprintf(stderr, "%d: ghost_comm send to %d (%d bytes)\n",
                            this_node, node, n_s_buffer));
        if (n_s_buffer > 0) {
          gc.mpi_comm.send(node, REQ_GHOST_SEND, s_buffer, n_s_buffer);
        }

        assert((data_parts & GHOSTTRANS_PROPRTS) or (s_bondbuffer.empty()));

        if (not s_bondbuffer.empty()) {
          gc.mpi_comm.send(node, REQ_GHOST_SEND, s_bondbuffer.data(),
                           s_bondbuffer.size());
        }
        break;
      }
      case GHOST_BCST:
        GHOST_TRACE(fprintf(stderr, "%d: ghost_comm bcast from %d (%d bytes)\n",
                            this_node, node,
                            (node == this_node) ? n_s_buffer : n_r_buffer));
        if (node == this_node) {
          if (n_s_buffer > 0) {
            boost::mpi::broadcast(gc.mpi_comm, s_buffer, n_s_buffer, node);
          }

          assert((data_parts & GHOSTTRANS_PROPRTS) or (s_bondbuffer.empty()));
          if (not s_bondbuffer.empty()) {
            boost::mpi::broadcast(gc.mpi_comm, s_bondbuffer.data(),
                                  s_bondbuffer.size(), node);
          }
        } else {
          if (n_r_buffer > 0) {
            boost::mpi::broadcast(gc.mpi_comm, r_buffer, n_r_buffer, node);
          }
          if (data_parts & GHOSTTRANS_PROPRTS) {
            int n_bonds = *(int *)(r_buffer + n_r_buffer - sizeof(int));
            if (n_bonds) {
              r_bondbuffer.resize(n_bonds);
              boost::mpi::broadcast(gc.mpi_comm, r_bondbuffer.data(),
                                    r_bondbuffer.size(), node);
            }
          }
        }
        break;
      case GHOST_RDCE: {
        GHOST_TRACE(fprintf(stderr, "%d: ghost_comm reduce to %d (%d bytes)\n",
                            this_node, node, n_s_buffer));

        if (node == this_node) {
          if (n_s_buffer > 0) {

            MPI_Reduce(reinterpret_cast<double *>(s_buffer),
                       reinterpret_cast<double *>(r_buffer),
                       n_s_buffer / sizeof(double), MPI_DOUBLE, MPI_SUM, node,
                       gc.mpi_comm);
          }
        } else {
          if (n_s_buffer > 0) {
            MPI_Reduce(reinterpret_cast<double *>(s_buffer), nullptr,
                       n_s_buffer / sizeof(double), MPI_DOUBLE, MPI_SUM, node,
                       gc.mpi_comm);
          }
        }
      } break;
      }
      GHOST_TRACE(fprintf(stderr, "%d: ghost_comm done\n", this_node));

      /* recv op; write back data directly, if no PSTSTORE delay is requested.
       */
      if (is_recv_op(comm_type, node)) {
        if (!poststore) {
          /* forces have to be added, the rest overwritten. Exception is RDCE,
             where the addition is integrated into the communication. */
          if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
            add_forces_from_recv_buffer(gcn->part_lists, {reinterpret_cast<ParticleForce *>(r_buffer),
                                                          n_r_buffer / sizeof(ParticleForce)});
          else
            put_recv_buffer(gcn, data_parts);
        } else {
          GHOST_TRACE(fprintf(
              stderr, "%d: ghost_comm delaying operation %d, recv from %d\n",
              this_node, n, node));
        }
      } else {
        /* send op; write back delayed data from last recv, when this was a
         * prefetch send. */
        if (poststore) {
          /* find previous action where we recv and which has PSTSTORE set */
          for (auto jt = make_reverse_iterator(it); jt != gc.comm.rend();
               ++jt) {
            GhostCommunication *gcn2 = &(*jt);
            int comm_type2 = gcn2->type;
            int poststore2 = gcn2->poststore;
            int node2 = gcn2->node;
            if (is_recv_op(comm_type2, node2) && poststore2) {
              assert(n_r_buffer == calc_transmit_size(gcn2->part_lists, data_parts));
              /* as above */
              if (data_parts == GHOSTTRANS_FORCE && comm_type != GHOST_RDCE)
                add_forces_from_recv_buffer(gcn2->part_lists, {reinterpret_cast<ParticleForce *>(r_buffer),
                                                               n_r_buffer / sizeof(ParticleForce)});
              else
                put_recv_buffer(gcn2, data_parts);
              break;
            }
          }
        }
      }
    }
  }
}
