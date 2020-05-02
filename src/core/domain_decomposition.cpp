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
 *  Implementation of domain_decomposition.hpp.
 */

#include "domain_decomposition.hpp"
#include "errorhandling.hpp"
#include "ghosts.hpp"

#include <utils/index.hpp>
#include <utils/mpi/sendrecv.hpp>
using Utils::get_linear_index;
#include <utils/mpi/cart_comm.hpp>

#include "cells.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/range/algorithm/reverse.hpp>

/************************************************/
/** \name Variables */
/************************************************/
/*@{*/

DomainDecomposition dd;

/*@}*/

/************************************************************/
/** \name Private Functions */
/************************************************************/

namespace {
/* Calc the ghost shift vector for dim dir in direction lr */
Utils::Vector3d shift(BoxGeometry const &box, LocalBox<double> const &local_box,
                      int dir, int lr) {
  Utils::Vector3d ret{};

  /* Shift is non-zero only in periodic directions, if we are at the box
   * boundary */
  ret[dir] = box.periodic(dir) * local_box.boundary()[2 * dir + lr] *
             box.length()[dir];

  return ret;
}
} // namespace

/** Create communicators for cell structure domain decomposition. (see \ref
 *  GhostCommunicator)
 */
GhostCommunicator dd_prepare_comm(const BoxGeometry &box_geo,
                                  const LocalBox<double> &local_geo) {
  int dir, lr, i, cnt, n_comm_cells[3];
  Utils::Vector3i lc{}, hc{}, done{};

  auto const comm_info = Utils::Mpi::cart_get<3>(dd.comm);
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(dd.comm);

  /* calculate number of communications */
  size_t num = 0;
  for (dir = 0; dir < 3; dir++) {
    for (lr = 0; lr < 2; lr++) {
      /* No communication for border of non periodic direction */
      if (comm_info.dims[dir] == 1)
        num++;
      else
        num += 2;
    }
  }

  /* prepare communicator */
  auto ghost_comm = GhostCommunicator{dd.comm, num};

  /* number of cells to communicate in a direction */
  n_comm_cells[0] = dd.cell_grid[1] * dd.cell_grid[2];
  n_comm_cells[1] = dd.cell_grid[2] * dd.ghost_cell_grid[0];
  n_comm_cells[2] = dd.ghost_cell_grid[0] * dd.ghost_cell_grid[1];

  cnt = 0;
  /* direction loop: x, y, z */
  for (dir = 0; dir < 3; dir++) {
    lc[(dir + 1) % 3] = 1 - done[(dir + 1) % 3];
    lc[(dir + 2) % 3] = 1 - done[(dir + 2) % 3];
    hc[(dir + 1) % 3] = dd.cell_grid[(dir + 1) % 3] + done[(dir + 1) % 3];
    hc[(dir + 2) % 3] = dd.cell_grid[(dir + 2) % 3] + done[(dir + 2) % 3];
    /* lr loop: left right */
    /* here we could in principle build in a one sided ghost
       communication, simply by taking the lr loop only over one
       value */
    for (lr = 0; lr < 2; lr++) {
      if (comm_info.dims[dir] == 1) {
        /* just copy cells on a single node */
        ghost_comm.communications[cnt].type = GHOST_LOCL;
        ghost_comm.communications[cnt].node = dd.comm.rank();

        /* Buffer has to contain Send and Recv cells -> factor 2 */
        ghost_comm.communications[cnt].part_lists.resize(2 * n_comm_cells[dir]);
        /* prepare folding of ghost positions */
        ghost_comm.communications[cnt].shift =
            shift(box_geo, local_geo, dir, lr);

        /* fill send ghost_comm cells */
        lc[dir] = hc[dir] = 1 + lr * (dd.cell_grid[dir] - 1);

        dd.fill_comm_cell_lists(
            ghost_comm.communications[cnt].part_lists.data(), lc, hc);

        /* fill recv ghost_comm cells */
        lc[dir] = hc[dir] = 0 + (1 - lr) * (dd.cell_grid[dir] + 1);

        /* place receive cells after send cells */
        dd.fill_comm_cell_lists(
            &ghost_comm.communications[cnt].part_lists[n_comm_cells[dir]], lc,
            hc);

        cnt++;
      } else {
        /* i: send/recv loop */
        for (i = 0; i < 2; i++) {
          if ((comm_info.coords[dir] + i) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_SEND;
            ghost_comm.communications[cnt].node = node_neighbors[2 * dir + lr];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);
            /* prepare folding of ghost positions */
            ghost_comm.communications[cnt].shift =
                shift(box_geo, local_geo, dir, lr);

            lc[dir] = hc[dir] = 1 + lr * (dd.cell_grid[dir] - 1);

            dd.fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
          if ((comm_info.coords[dir] + (1 - i)) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_RECV;
            ghost_comm.communications[cnt].node =
                node_neighbors[2 * dir + (1 - lr)];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);

            lc[dir] = hc[dir] = (1 - lr) * (dd.cell_grid[dir] + 1);

            dd.fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
        }
      }
      done[dir] = 1;
    }
  }

  return ghost_comm;
}

namespace {
/** Revert the order of a communicator: After calling this the
 *  communicator is working in reverted order with exchanged
 *  communication types GHOST_SEND <-> GHOST_RECV.
 */
void revert_comm_order(GhostCommunicator &comm) {
  /* revert order */
  boost::reverse(comm.communications);

  /* exchange SEND/RECV */
  for (auto &c : comm.communications) {
    if (c.type == GHOST_SEND)
      c.type = GHOST_RECV;
    else if (c.type == GHOST_RECV)
      c.type = GHOST_SEND;
    else if (c.type == GHOST_LOCL) {
      boost::reverse(c.part_lists);
    }
  }
}

/** Of every two communication rounds, set the first receivers to prefetch and
 *  poststore
 */
void assign_prefetches(GhostCommunicator &comm) {
  for (auto it = comm.communications.begin(); it != comm.communications.end();
       it += 2) {
    auto next = std::next(it);
    if (it->type == GHOST_RECV && next->type == GHOST_SEND) {
      it->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
      next->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
    }
  }
}
} // namespace

/** update the 'shift' member of those GhostCommunicators, which use
 *  that value to speed up the folding process of its ghost members
 *  (see \ref dd_prepare_comm for the original).
 */
void dd_update_communicators_w_boxl() {
  auto const cart_info = Utils::Mpi::cart_get<3>(dd.comm);

  int cnt = 0;
  /* direction loop: x, y, z */
  for (int dir = 0; dir < 3; dir++) {
    /* lr loop: left right */
    for (int lr = 0; lr < 2; lr++) {
      if (cart_info.dims[dir] == 1) {
        /* prepare folding of ghost positions */
        if (dd.local_geo.boundary()[2 * dir + lr] != 0) {
          dd.m_exchange_ghosts_comm.communications[cnt].shift =
              shift(dd.box_geo, dd.local_geo, dir, lr);
        }
        cnt++;
      } else {
        /* i: send/recv loop */
        for (int i = 0; i < 2; i++) {
          if ((cart_info.coords[dir] + i) % 2 == 0) {
            /* prepare folding of ghost positions */
            if (dd.local_geo.boundary()[2 * dir + lr] != 0) {
              dd.m_exchange_ghosts_comm.communications[cnt].shift =
                  shift(dd.box_geo, dd.local_geo, dir, lr);
            }
            cnt++;
          }
          if ((cart_info.coords[dir] + (1 - i)) % 2 == 0) {
            cnt++;
          }
        }
      }
    }
  }
}

/*************************************************/

/*@}*/

/************************************************************/
/* Public Functions */
/************************************************************/

bool dd_on_geometry_change(bool fast, double range, const BoxGeometry &box_geo,
                           const LocalBox<double> &local_geo) {
  /* check that the CPU domains are still sufficiently large. */
  for (int i = 0; i < 3; i++)
    if (local_geo.length()[i] < range) {
      runtimeErrorMsg() << "local box length " << local_geo.length()[i]
                        << " in direction " << i
                        << " is smaller than"
                           "interaction radius "
                        << range;
      return false;
    }

  double min_cell_size =
      std::min(std::min(dd.cell_size[0], dd.cell_size[1]), dd.cell_size[2]);

  /* If new box length leads to too small cells, redo cell structure
   using smaller number of cells. If we are not in a hurry, check if we can
   maybe optimize the cell system by using smaller cells. */
  auto const re_init = (range > min_cell_size) or ((not fast) and [&]() {
                         for (int i = 0; i < 3; i++) {
                           auto const poss_size = static_cast<int>(
                               floor(local_geo.length()[i] / range));
                           if (poss_size > dd.cell_grid[i])
                             return true;
                         }
                         return false;
                       }());

  if (re_init) {
    return false;
  }

  dd.box_geo = box_geo;
  dd.local_geo = local_geo;

  /* otherwise, re-set our geometrical dimensions which have changed */
  for (int i = 0; i < 3; i++) {
    dd.cell_size[i] = local_geo.length()[i] / (double)dd.cell_grid[i];
    dd.inv_cell_size[i] = 1.0 / dd.cell_size[i];
  }

  dd_update_communicators_w_boxl();

  return true;
}

/************************************************************/
void dd_topology_init(const boost::mpi::communicator &comm, double range,
                      const BoxGeometry &box_geo,
                      const LocalBox<double> &local_geo) {
  dd.comm = comm;
  auto const cart_info = Utils::Mpi::cart_get<3>(dd.comm);

  for (int i = 0; i < 3; i++) {
    if (dd.fully_connected[i] and cart_info.dims[i] != 1) {
      runtimeErrorMsg()
          << "Node grid not compatible with fully_connected property";
    }
  }

  dd.box_geo = box_geo;
  dd.local_geo = local_geo;

  /* set up new domain decomposition cell structure */
  dd.create_cell_grid(range);

  /* setup cell neighbors */
  dd.init_cell_interactions();

  /* mark local and ghost cells */
  dd.mark_cells();

  /* create communicators */
  dd.m_exchange_ghosts_comm = dd_prepare_comm(box_geo, local_geo);
  dd.m_collect_ghost_force_comm = dd_prepare_comm(box_geo, local_geo);

  /* collect forces has to be done in reverted order! */
  revert_comm_order(dd.m_collect_ghost_force_comm);

  assign_prefetches(dd.m_exchange_ghosts_comm);
  assign_prefetches(dd.m_collect_ghost_force_comm);

  cell_structure.type = CELL_STRUCTURE_DOMDEC;
  cell_structure.m_decomposition = std::addressof(dd);
}

void dd_exchange_and_sort_particles(int global, ParticleList *pl,
                                    std::vector<Cell *> &modified_cells) {
  dd.resort(global, (assert(pl), *pl), modified_cells);
}

/************************************************************/
