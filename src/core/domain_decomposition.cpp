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

/** Returns pointer to the cell which corresponds to the position if the
 *  position is in the nodes spatial domain otherwise a nullptr pointer.
 */
Cell *dd_save_position_to_cell(const Utils::Vector3d &pos);

/************************************************/
/** \name Variables */
/************************************************/
/*@{*/

DomainDecomposition dd;

/** Maximal number of cells per node. In order to avoid memory
 *  problems due to the cell grid one has to specify the maximal
 *  number of \ref cells::cells. If the number of cells is larger
 *  than max_num_cells the cell grid is reduced.
 *  max_num_cells has to be larger than 27, e.g. one inner cell.
 */
constexpr int max_num_cells = 32768;

/*@}*/

/************************************************************/
/** \name Private Functions */
/************************************************************/
/*@{*/

/**
 *  @brief Calculate cell grid dimensions, cell sizes and number of cells.
 *
 *  Calculates the cell grid, based on \ref local_geo and \p range.
 *  If the number of cells is larger than \ref max_num_cells,
 *  it increases max_range until the number of cells is
 *  smaller or equal \ref max_num_cells. It sets:
 *  \ref DomainDecomposition::cell_grid,
 *  \ref DomainDecomposition::ghost_cell_grid,
 *  \ref DomainDecomposition::cell_size, and
 *  \ref DomainDecomposition::inv_cell_size.
 *
 *  @param range Required interacting range. All pairs closer
 *         than this distance are found.
 */
void dd_create_cell_grid(double range) {
  auto const cart_info = Utils::Mpi::cart_get<3>(dd.comm);

  int i, n_local_cells, new_cells;
  double cell_range[3];

  /* initialize */
  cell_range[0] = cell_range[1] = cell_range[2] = range;

  /* Min num cells can not be smaller than calc_processor_min_num_cells. */
  int min_num_cells = calc_processor_min_num_cells(cart_info.dims);

  if (range <= 0.) {
    /* this is the non-interacting case */
    const int cells_per_dir = std::ceil(std::pow(min_num_cells, 1. / 3.));

    dd.cell_grid[0] = cells_per_dir;
    dd.cell_grid[1] = cells_per_dir;
    dd.cell_grid[2] = cells_per_dir;

    n_local_cells = dd.cell_grid[0] * dd.cell_grid[1] * dd.cell_grid[2];
  } else {
    /* Calculate initial cell grid */
    double volume = dd.local_geo.length()[0];
    for (i = 1; i < 3; i++)
      volume *= dd.local_geo.length()[i];
    double scale = pow(max_num_cells / volume, 1. / 3.);
    for (i = 0; i < 3; i++) {
      /* this is at least 1 */
      dd.cell_grid[i] = (int)ceil(dd.local_geo.length()[i] * scale);
      cell_range[i] = dd.local_geo.length()[i] / dd.cell_grid[i];

      if (cell_range[i] < range) {
        /* ok, too many cells for this direction, set to minimum */
        dd.cell_grid[i] = (int)floor(dd.local_geo.length()[i] / range);
        if (dd.cell_grid[i] < 1) {
          runtimeErrorMsg() << "interaction range " << range << " in direction "
                            << i << " is larger than the local box size "
                            << dd.local_geo.length()[i];
          dd.cell_grid[i] = 1;
        }
        cell_range[i] = dd.local_geo.length()[i] / dd.cell_grid[i];
      }
    }

    /* It may be necessary to asymmetrically assign the scaling to the
       coordinates, which the above approach will not do.
       For a symmetric box, it gives a symmetric result. Here we correct that.
       */
    for (;;) {
      n_local_cells = dd.cell_grid[0] * dd.cell_grid[1] * dd.cell_grid[2];

      /* done */
      if (n_local_cells <= max_num_cells)
        break;

      /* find coordinate with the smallest cell range */
      int min_ind = 0;
      double min_size = cell_range[0];

      for (i = 1; i < 3; i++) {
        if (dd.cell_grid[i] > 1 && cell_range[i] < min_size) {
          min_ind = i;
          min_size = cell_range[i];
        }
      }

      dd.cell_grid[min_ind]--;
      cell_range[min_ind] =
          dd.local_geo.length()[min_ind] / dd.cell_grid[min_ind];
    }

    /* sanity check */
    if (n_local_cells < min_num_cells) {
      runtimeErrorMsg()
          << "number of cells " << n_local_cells << " is smaller than minimum "
          << min_num_cells
          << " (interaction range too large or min_num_cells too large)";
    }
  }

  /* quit program if unsuccessful */
  if (n_local_cells > max_num_cells) {
    runtimeErrorMsg() << "no suitable cell grid found ";
  }

  auto const node_pos = cart_info.coords;

  /* now set all dependent variables */
  new_cells = 1;
  for (i = 0; i < 3; i++) {
    dd.ghost_cell_grid[i] = dd.cell_grid[i] + 2;
    new_cells *= dd.ghost_cell_grid[i];
    dd.cell_size[i] = dd.local_geo.length()[i] / (double)dd.cell_grid[i];
    dd.inv_cell_size[i] = 1.0 / dd.cell_size[i];
    dd.cell_offset[i] = node_pos[i] * dd.cell_grid[i];
  }
  cell_structure.max_range = dd.cell_size;

  /* allocate cell array and cell pointer arrays */
  cells.resize(new_cells);
  cell_structure.m_local_cells.resize(n_local_cells);
  cell_structure.m_ghost_cells.resize(new_cells - n_local_cells);
}

/** Fill local_cells list and ghost_cells list for use with domain
 *  decomposition.  \ref cells::cells is assumed to be a 3d grid with size
 *  \ref DomainDecomposition::ghost_cell_grid.
 */
void dd_mark_cells() {
  int cnt_c = 0, cnt_l = 0, cnt_g = 0;

  for (int o = 0; o < dd.ghost_cell_grid[2]; o++)
    for (int n = 0; n < dd.ghost_cell_grid[1]; n++)
      for (int m = 0; m < dd.ghost_cell_grid[0]; m++) {
        if ((m > 0 && m < dd.ghost_cell_grid[0] - 1 && n > 0 &&
             n < dd.ghost_cell_grid[1] - 1 && o > 0 &&
             o < dd.ghost_cell_grid[2] - 1))
          cell_structure.m_local_cells[cnt_l++] = &cells[cnt_c++];
        else
          cell_structure.m_ghost_cells[cnt_g++] = &cells[cnt_c++];
      }
}

/** Fill a communication cell pointer list. Fill the cell pointers of
 *  all cells which are inside a rectangular subgrid of the 3D cell
 *  grid (\ref DomainDecomposition::ghost_cell_grid) starting from the
 *  lower left corner lc up to the high top corner hc. The cell
 *  pointer list part_lists must already be large enough.
 *  \param part_lists  List of cell pointers to store the result.
 *  \param lc          lower left corner of the subgrid.
 *  \param hc          high up corner of the subgrid.
 */
int dd_fill_comm_cell_lists(ParticleList **part_lists,
                            Utils::Vector3i const &lc,
                            Utils::Vector3i const &hc) {
  /* sanity check */
  for (int i = 0; i < 3; i++) {
    if (lc[i] < 0 || lc[i] >= dd.ghost_cell_grid[i])
      return 0;
    if (hc[i] < 0 || hc[i] >= dd.ghost_cell_grid[i])
      return 0;
    if (lc[i] > hc[i])
      return 0;
  }

  int c = 0;
  for (int o = lc[0]; o <= hc[0]; o++)
    for (int n = lc[1]; n <= hc[1]; n++)
      for (int m = lc[2]; m <= hc[2]; m++) {
        auto const i =
            get_linear_index(o, n, m,
                             {dd.ghost_cell_grid[0], dd.ghost_cell_grid[1],
                              dd.ghost_cell_grid[2]});

        part_lists[c] = &(cells[i].particles());
        c++;
      }
  return c;
}

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

        dd_fill_comm_cell_lists(
            ghost_comm.communications[cnt].part_lists.data(), lc, hc);

        /* fill recv ghost_comm cells */
        lc[dir] = hc[dir] = 0 + (1 - lr) * (dd.cell_grid[dir] + 1);

        /* place receive cells after send cells */
        dd_fill_comm_cell_lists(
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

            dd_fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
          if ((comm_info.coords[dir] + (1 - i)) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_RECV;
            ghost_comm.communications[cnt].node =
                node_neighbors[2 * dir + (1 - lr)];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);

            lc[dir] = hc[dir] = (1 - lr) * (dd.cell_grid[dir] + 1);

            dd_fill_comm_cell_lists(
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

/** Revert the order of a communicator: After calling this the
 *  communicator is working in reverted order with exchanged
 *  communication types GHOST_SEND <-> GHOST_RECV.
 */
void dd_revert_comm_order(GhostCommunicator *comm) {
  /* revert order */
  boost::reverse(comm->communications);

  /* exchange SEND/RECV */
  for (auto &c : comm->communications) {
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
void dd_assign_prefetches(GhostCommunicator *comm) {
  for (auto it = comm->communications.begin(); it != comm->communications.end();
       it += 2) {
    auto next = std::next(it);
    if (it->type == GHOST_RECV && next->type == GHOST_SEND) {
      it->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
      next->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
    }
  }
}

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
          cell_structure.exchange_ghosts_comm.communications[cnt].shift =
              shift(dd.box_geo, dd.local_geo, dir, lr);
        }
        cnt++;
      } else {
        /* i: send/recv loop */
        for (int i = 0; i < 2; i++) {
          if ((cart_info.coords[dir] + i) % 2 == 0) {
            /* prepare folding of ghost positions */
            if (dd.local_geo.boundary()[2 * dir + lr] != 0) {
              cell_structure.exchange_ghosts_comm.communications[cnt].shift =
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

/** Init cell interactions for cell system domain decomposition.
 * initializes the interacting neighbor cell list of a cell The
 * created list of interacting neighbor cells is used by the Verlet
 * algorithm (see verlet.cpp) to build the verlet lists.
 */
void dd_init_cell_interactions(const Utils::Vector3i &grid) {
  int m, n, o, p, q, r, ind1, ind2;

  for (int i = 0; i < 3; i++) {
    if (dd.fully_connected[i] and grid[i] != 1) {
      runtimeErrorMsg()
          << "Node grid not compatible with fully_connected property";
    }
  }

  /* loop all local cells */
  for (o = 1; o < dd.cell_grid[2] + 1; o++)
    for (n = 1; n < dd.cell_grid[1] + 1; n++)
      for (m = 1; m < dd.cell_grid[0] + 1; m++) {

        ind1 = get_linear_index(m, n, o,
                                {dd.ghost_cell_grid[0], dd.ghost_cell_grid[1],
                                 dd.ghost_cell_grid[2]});

        std::vector<Cell *> red_neighbors;
        std::vector<Cell *> black_neighbors;

        /* loop all neighbor cells */
        int lower_index[3] = {m - 1, n - 1, o - 1};
        int upper_index[3] = {m + 1, n + 1, o + 1};

        for (int i = 0; i < 3; i++) {
          if (dd.fully_connected[i]) {
            lower_index[i] = 0;
            upper_index[i] = dd.ghost_cell_grid[i] - 1;
          }
        }

        for (p = lower_index[2]; p <= upper_index[2]; p++)
          for (q = lower_index[1]; q <= upper_index[1]; q++)
            for (r = lower_index[0]; r <= upper_index[0]; r++) {
              ind2 = get_linear_index(r, q, p,
                                      {dd.ghost_cell_grid[0],
                                       dd.ghost_cell_grid[1],
                                       dd.ghost_cell_grid[2]});
              if (ind2 > ind1) {
                red_neighbors.push_back(&cells[ind2]);
              } else {
                black_neighbors.push_back(&cells[ind2]);
              }
            }
        cells[ind1].m_neighbors =
            Neighbors<Cell *>(red_neighbors, black_neighbors);
      }
}

/*************************************************/

/** Returns pointer to the cell which corresponds to the position if the
 *  position is in the nodes spatial domain otherwise a nullptr pointer.
 */
Cell *dd_save_position_to_cell(const Utils::Vector3d &pos) {
  int cpos[3];

  for (int i = 0; i < 3; i++) {
    cpos[i] = static_cast<int>(std::floor(pos[i] * dd.inv_cell_size[i])) + 1 -
              dd.cell_offset[i];

    /* particles outside our box. Still take them if
       nonperiodic boundary. We also accept the particle if we are at
       the box boundary, and the particle is within the box. In this case
       the particle belongs here and could otherwise potentially be dismissed
       due to rounding errors. */
    if (cpos[i] < 1) {
      if ((!dd.box_geo.periodic(i) or (pos[i] >= dd.box_geo.length()[i])) &&
          dd.local_geo.boundary()[2 * i])
        cpos[i] = 1;
      else
        return nullptr;
    } else if (cpos[i] > dd.cell_grid[i]) {
      if ((!dd.box_geo.periodic(i) or (pos[i] < dd.box_geo.length()[i])) &&
          dd.local_geo.boundary()[2 * i + 1])
        cpos[i] = dd.cell_grid[i];
      else
        return nullptr;
    }
  }

  auto const ind = get_linear_index(
      cpos[0], cpos[1], cpos[2],
      {dd.ghost_cell_grid[0], dd.ghost_cell_grid[1], dd.ghost_cell_grid[2]});
  return &(cells[ind]);
}

/*@}*/

/************************************************************/
/* Public Functions */
/************************************************************/

void dd_on_geometry_change(bool fast, double range, const BoxGeometry &box_geo,
                           const LocalBox<double> &local_geo) {
  /* check that the CPU domains are still sufficiently large. */
  for (int i = 0; i < 3; i++)
    if (local_geo.length()[i] < range) {
      runtimeErrorMsg() << "local box length " << local_geo.length()[i]
                        << " in direction " << i
                        << " is smaller than"
                           "interaction radius "
                        << range;
    }

  dd.box_geo = box_geo;
  dd.local_geo = local_geo;

  /* otherwise, re-set our geometrical dimensions which have changed */
  for (int i = 0; i < 3; i++) {
    dd.cell_size[i] = local_geo.length()[i] / (double)dd.cell_grid[i];
    dd.inv_cell_size[i] = 1.0 / dd.cell_size[i];
  }

  double min_cell_size =
      std::min(std::min(dd.cell_size[0], dd.cell_size[1]), dd.cell_size[2]);

  if (range > min_cell_size) {
    /* if new box length leads to too small cells, redo cell structure
       using smaller number of cells. */
    cells_re_init(CELL_STRUCTURE_DOMDEC, range);
    return;
  }

  /* If we are not in a hurry, check if we can maybe optimize the cell
     system by using smaller cells. */
  if (not fast) {
    int i;
    for (i = 0; i < 3; i++) {
      auto poss_size = (int)floor(local_geo.length()[i] / range);
      if (poss_size > dd.cell_grid[i])
        break;
    }
    if (i < 3) {
      /* new range/box length allow smaller cells, redo cell structure,
         possibly using smaller number of cells. */
      cells_re_init(CELL_STRUCTURE_DOMDEC, range);
      return;
    }
  }
  dd_update_communicators_w_boxl();
}

int calc_processor_min_num_cells(const Utils::Vector3i &grid) {
  int i, min = 1;
  /* the minimal number of cells can be lower if there are at least two nodes
     serving a direction,
     since this also ensures that the cell size is at most half the box length.
     However, if there is
     only one processor for a direction, there have to be at least two cells for
     this direction. */
  for (i = 0; i < 3; i++)
    if (grid[i] == 1)
      min *= 2;
  return min;
}

/************************************************************/
void dd_topology_init(const boost::mpi::communicator &comm, double range,
                      const BoxGeometry &box_geo,
                      const LocalBox<double> &local_geo) {
  dd.comm = comm;
  auto const cart_info = Utils::Mpi::cart_get<3>(dd.comm);

  dd.box_geo = box_geo;
  dd.local_geo = local_geo;

  cell_structure.type = CELL_STRUCTURE_DOMDEC;
  cell_structure.particle_to_cell = [](const Particle &p) {
    return dd_save_position_to_cell(p.r.p);
  };

  /* set up new domain decomposition cell structure */
  dd_create_cell_grid(range);
  /* mark cells */
  dd_mark_cells();

  /* create communicators */
  cell_structure.exchange_ghosts_comm = dd_prepare_comm(box_geo, local_geo);
  cell_structure.collect_ghost_force_comm = dd_prepare_comm(box_geo, local_geo);

  /* collect forces has to be done in reverted order! */
  dd_revert_comm_order(&cell_structure.collect_ghost_force_comm);

  dd_assign_prefetches(&cell_structure.exchange_ghosts_comm);
  dd_assign_prefetches(&cell_structure.collect_ghost_force_comm);

  dd_init_cell_interactions(cart_info.dims);
}

namespace {

/**
 * @brief Move particles into the cell system if it belongs to this node.
 *
 * Moves all particles from src into the local cell
 * system if they do belong here. Otherwise the
 * particles are moved into rest.
 *
 * @param src Particles to move.
 * @param rest Output list for left-over particles.
 * @param modified_cells Local cells that were touched.
 */
void move_if_local(ParticleList &src, ParticleList &rest,
                   std::vector<Cell *> &modified_cells) {
  for (auto &part : src) {
    assert(cell_structure.get_local_particle(part.p.identity) == nullptr);

    auto target_cell = dd_save_position_to_cell(part.r.p);

    if (target_cell) {
      target_cell->particles().insert(std::move(part));
      modified_cells.push_back(target_cell);
    } else {
      rest.insert(std::move(part));
    }
  }

  src.clear();
}

/**
 * @brief Split particle list by direction.
 *
 * Moves all particles from src into left
 * and right depending if they belong to
 * the left or right side from local node
 * in direction dir.
 *
 * @param src Particles to sort.
 * @param left Particles that should go to the left
 * @param right Particles that should go to the right
 * @param dir Direction to consider.
 */
void move_left_or_right(ParticleList &src, ParticleList &left,
                        ParticleList &right, int dir) {
  for (auto it = src.begin(); it != src.end();) {
    assert(cell_structure.get_local_particle(it->p.identity) == nullptr);

    if ((get_mi_coord(it->r.p[dir], dd.local_geo.my_left()[dir],
                      dd.box_geo.length()[dir],
                      dd.box_geo.periodic(dir)) < 0.0) and
        (dd.box_geo.periodic(dir) || (dd.local_geo.boundary()[2 * dir] == 0))) {
      left.insert(std::move(*it));
      it = src.erase(it);
    } else if ((get_mi_coord(it->r.p[dir], dd.local_geo.my_right()[dir],
                             dd.box_geo.length()[dir],
                             dd.box_geo.periodic(dir)) >= 0.0) and
               (dd.box_geo.periodic(dir) ||
                (dd.local_geo.boundary()[2 * dir + 1] == 0))) {
      right.insert(std::move(*it));
      it = src.erase(it);
    } else {
      ++it;
    }
  }
} // namespace

void exchange_neighbors(ParticleList *pl, std::vector<Cell *> &modified_cells) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(dd.comm);
  static ParticleList send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;

  for (int dir = 0; dir < 3; dir++) {
    /* Single node direction, no action needed. */
    if (Utils::Mpi::cart_get<3>(dd.comm).dims[dir] == 1) {
      continue;
      /* In this (common) case left and right neighbors are
         the same, and we need only one communication */
    }
    if (Utils::Mpi::cart_get<3>(dd.comm).dims[dir] == 2) {
      move_left_or_right(*pl, send_buf_l, send_buf_l, dir);

      Utils::Mpi::sendrecv(dd.comm, node_neighbors[2 * dir], 0, send_buf_l,
                           node_neighbors[2 * dir], 0, recv_buf_l);

      send_buf_l.clear();
    } else {
      using boost::mpi::request;
      using Utils::Mpi::isendrecv;

      move_left_or_right(*pl, send_buf_l, send_buf_r, dir);

      auto req_l = isendrecv(dd.comm, node_neighbors[2 * dir], 0, send_buf_l,
                             node_neighbors[2 * dir], 0, recv_buf_l);
      auto req_r =
          isendrecv(dd.comm, node_neighbors[2 * dir + 1], 0, send_buf_r,
                    node_neighbors[2 * dir + 1], 0, recv_buf_r);

      std::array<request, 4> reqs{{req_l[0], req_l[1], req_r[0], req_r[1]}};
      boost::mpi::wait_all(reqs.begin(), reqs.end());

      send_buf_l.clear();
      send_buf_r.clear();
    }

    move_if_local(recv_buf_l, *pl, modified_cells);
    move_if_local(recv_buf_r, *pl, modified_cells);
  }
}
} // namespace

void dd_exchange_and_sort_particles(int global, ParticleList *pl,
                                    std::vector<Cell *> &modified_cells) {
  if (global) {
    auto const grid = Utils::Mpi::cart_get<3>(dd.comm).dims;

    /* Worst case we need grid - 1 rounds per direction.
     * This correctly implies that if there is only one node,
     * no action should be taken. */
    int rounds_left = grid[0] + grid[1] + grid[2] - 3;
    for (; rounds_left > 0; rounds_left--) {
      exchange_neighbors(pl, modified_cells);

      auto left_over =
          boost::mpi::all_reduce(dd.comm, pl->size(), std::plus<size_t>());

      if (left_over == 0) {
        break;
      }
    }
  } else {
    exchange_neighbors(pl, modified_cells);
  }
}

/************************************************************/
