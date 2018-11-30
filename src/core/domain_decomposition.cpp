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
 *
 *  This file contains everything related to the cell system: domain
 * decomposition.
 *  See also \ref domain_decomposition.hpp
 */

#include "domain_decomposition.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include "utils/mpi/sendrecv.hpp"
#include "utils/serialization/ParticleList.hpp"

#include "initialize.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/range/algorithm.hpp>

/** Returns pointer to the cell which corresponds to the position if
    the position is in the nodes spatial domain otherwise a nullptr
    pointer. */
Cell *dd_save_position_to_cell(const Vector3d &pos);

/************************************************/
/** \name Defines */
/************************************************/
/*@{*/

/** half the number of cell neighbors in 3 Dimensions. */
#define CELLS_MAX_NEIGHBORS 14

/*@}*/

/************************************************/
/** \name Variables */
/************************************************/
/*@{*/

DomainDecomposition dd;

int max_num_cells = CELLS_MAX_NUM_CELLS;
int min_num_cells = 1;
double max_skin = 0.0;

/*@}*/

/************************************************************/
/** \name Private Functions */
/************************************************************/
/*@{*/

/** Convenient replace for loops over all cells. */
#define DD_CELLS_LOOP(m, n, o)                                                 \
  for (o = 0; o < dd.ghost_cell_grid[2]; o++)                                  \
    for (n = 0; n < dd.ghost_cell_grid[1]; n++)                                \
      for (m = 0; m < dd.ghost_cell_grid[0]; m++)

/** Convenient replace for loops over Local cells. */
#define DD_LOCAL_CELLS_LOOP(m, n, o)                                           \
  for (o = 1; o < dd.cell_grid[2] + 1; o++)                                    \
    for (n = 1; n < dd.cell_grid[1] + 1; n++)                                  \
      for (m = 1; m < dd.cell_grid[0] + 1; m++)

/** Convenient replace for inner cell check. usage: if(DD_IS_LOCAL_CELL(m,n,o))
 * {...} */
#define DD_IS_LOCAL_CELL(m, n, o)                                              \
  (m > 0 && m < dd.ghost_cell_grid[0] - 1 && n > 0 &&                          \
   n < dd.ghost_cell_grid[1] - 1 && o > 0 && o < dd.ghost_cell_grid[2] - 1)

/** Convenient replace for ghost cell check. usage: if(DD_IS_GHOST_CELL(m,n,o))
 * {...} */
#define DD_IS_GHOST_CELL(m, n, o)                                              \
  (m == 0 || m == dd.ghost_cell_grid[0] - 1 || n == 0 ||                       \
   n >= dd.ghost_cell_grid[1] - 1 || o == 0 || o == dd.ghost_cell_grid[2] - 1)

/** Calculate cell grid dimensions, cell sizes and number of cells.
 *  Calculates the cell grid, based on \ref local_box_l and \ref
 *  max_range. If the number of cells is larger than \ref
 *  max_num_cells, it increases max_range until the number of cells is
 *  smaller or equal \ref max_num_cells. It sets: \ref
 *  DomainDecomposition::cell_grid, \ref
 *  DomainDecomposition::ghost_cell_grid, \ref
 *  DomainDecomposition::cell_size, \ref
 *  DomainDecomposition::inv_cell_size, and \ref n_cells.
 */
void dd_create_cell_grid() {
  int i, n_local_cells, new_cells;
  double cell_range[3];
  CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid: max_range %f\n",
                     this_node, max_range));
  CELL_TRACE(fprintf(
      stderr, "%d: dd_create_cell_grid: local_box %f-%f, %f-%f, %f-%f,\n",
      this_node, my_left[0], my_right[0], my_left[1], my_right[1], my_left[2],
      my_right[2]));

  /* initialize */
  cell_range[0] = cell_range[1] = cell_range[2] = max_range;

  if (max_range < ROUND_ERROR_PREC * box_l[0]) {
    /* this is the non-interacting case */
    const int cells_per_dir = std::ceil(std::pow(min_num_cells, 1. / 3.));

    dd.cell_grid[0] = cells_per_dir;
    dd.cell_grid[1] = cells_per_dir;
    dd.cell_grid[2] = cells_per_dir;

    n_local_cells = dd.cell_grid[0] * dd.cell_grid[1] * dd.cell_grid[2];
  } else {
    /* Calculate initial cell grid */
    double volume = local_box_l[0];
    for (i = 1; i < 3; i++)
      volume *= local_box_l[i];
    double scale = pow(max_num_cells / volume, 1. / 3.);
    for (i = 0; i < 3; i++) {
      /* this is at least 1 */
      dd.cell_grid[i] = (int)ceil(local_box_l[i] * scale);
      cell_range[i] = local_box_l[i] / dd.cell_grid[i];

      if (cell_range[i] < max_range) {
        /* ok, too many cells for this direction, set to minimum */
        dd.cell_grid[i] = (int)floor(local_box_l[i] / max_range);
        if (dd.cell_grid[i] < 1) {
          runtimeErrorMsg()
              << "interaction range " << max_range << " in direction " << i
              << " is larger than the local box size " << local_box_l[i];
          dd.cell_grid[i] = 1;
        }
        cell_range[i] = local_box_l[i] / dd.cell_grid[i];
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
      CELL_TRACE(fprintf(stderr,
                         "%d: minimal coordinate %d, size %f, grid %d\n",
                         this_node, min_ind, min_size, dd.cell_grid[min_ind]));

      dd.cell_grid[min_ind]--;
      cell_range[min_ind] = local_box_l[min_ind] / dd.cell_grid[min_ind];
    }
    CELL_TRACE(fprintf(stderr, "%d: final %d %d %d\n", this_node,
                       dd.cell_grid[0], dd.cell_grid[1], dd.cell_grid[2]));

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

  /* now set all dependent variables */
  new_cells = 1;
  for (i = 0; i < 3; i++) {
    dd.ghost_cell_grid[i] = dd.cell_grid[i] + 2;
    new_cells *= dd.ghost_cell_grid[i];
    dd.cell_size[i] = local_box_l[i] / (double)dd.cell_grid[i];
    dd.inv_cell_size[i] = 1.0 / dd.cell_size[i];
  }
  max_skin =
      std::min(std::min(dd.cell_size[0], dd.cell_size[1]), dd.cell_size[2]) -
      max_cut;

  /* allocate cell array and cell pointer arrays */
  realloc_cells(new_cells);
  realloc_cellplist(&local_cells, local_cells.n = n_local_cells);
  realloc_cellplist(&ghost_cells, ghost_cells.n = new_cells - n_local_cells);

  CELL_TRACE(fprintf(stderr,
                     "%d: dd_create_cell_grid, n_cells=%lu, local_cells.n=%d, "
                     "ghost_cells.n=%d, dd.ghost_cell_grid=(%d,%d,%d)\n",
                     this_node, (unsigned long)cells.size(), local_cells.n,
                     ghost_cells.n, dd.ghost_cell_grid[0],
                     dd.ghost_cell_grid[1], dd.ghost_cell_grid[2]));
}

/** Fill local_cells list and ghost_cells list for use with domain
    decomposition.  \ref cells::cells is assumed to be a 3d grid with size
    \ref DomainDecomposition::ghost_cell_grid . */
void dd_mark_cells() {
  int m, n, o, cnt_c = 0, cnt_l = 0, cnt_g = 0;

  DD_CELLS_LOOP(m, n, o) {

    if (DD_IS_LOCAL_CELL(m, n, o))
      local_cells.cell[cnt_l++] = &cells[cnt_c++];
    else
      ghost_cells.cell[cnt_g++] = &cells[cnt_c++];
  }
}

/** Fill a communication cell pointer list. Fill the cell pointers of
    all cells which are inside a rectangular subgrid of the 3D cell
    grid (\ref DomainDecomposition::ghost_cell_grid) starting from the
    lower left corner lc up to the high top corner hc. The cell
    pointer list part_lists must already be large enough.
    \param part_lists  List of cell pointers to store the result.
    \param lc          lower left corner of the subgrid.
    \param hc          high up corner of the subgrid.
 */
int dd_fill_comm_cell_lists(Cell **part_lists, int lc[3], int hc[3]) {
  int i, m, n, o, c = 0;
  /* sanity check */
  for (i = 0; i < 3; i++) {
    if (lc[i] < 0 || lc[i] >= dd.ghost_cell_grid[i])
      return 0;
    if (hc[i] < 0 || hc[i] >= dd.ghost_cell_grid[i])
      return 0;
    if (lc[i] > hc[i])
      return 0;
  }

  for (o = lc[0]; o <= hc[0]; o++)
    for (n = lc[1]; n <= hc[1]; n++)
      for (m = lc[2]; m <= hc[2]; m++) {
        i = get_linear_index(o, n, m,
                             {dd.ghost_cell_grid[0], dd.ghost_cell_grid[1],
                              dd.ghost_cell_grid[2]});
        CELL_TRACE(fprintf(stderr, "%d: dd_fill_comm_cell_list: add cell %d\n",
                           this_node, i));
        part_lists[c] = &cells[i];
        c++;
      }
  return c;
}

namespace {
    /**
     * @brief Helper to update a single component in an optional Vector3d.
     *
     * Updates a single component of an optional value, leaving the other
     * components unchanged if the optional had already a value, otherwise
     * initializing them to the default value.
     *
     * @param o the optional to update
     * @param value The updated value
     * @param i which one
     */
     void update_component(boost::optional<Vector3d> &o, double val, int i) {
         auto v = o.get_value_or({});
         v[i] = val;

        o = v;
     }
}

/** Create communicators for cell structure domain decomposition. (see \ref
 * GhostCommunicator) */
void dd_prepare_comm(GhostCommunicator *comm, int data_parts) {
  int i, num;
  int lc[3], hc[3], done[3] = {0, 0, 0};

  /* calculate number of communications */
  num = 0;
  for (int dir = 0; dir < 3; dir++) {
    for (int lr = 0; lr < 2; lr++) {
      /* No communication for border of non periodic direction */
      if (PERIODIC(dir) || (boundary[2 * dir + lr] == 0)) {
        if (node_grid[dir] == 1)
          num++;
        else
          num += 2;
      }
    }
  }

  /* prepare communicator */
  CELL_TRACE(fprintf(stderr,
                     "%d Create Communicator: prep_comm data_parts %d num %d\n",
                     this_node, data_parts, num));
  *comm = GhostCommunicator(data_parts, num);

  /* number of cells to communicate in a direction */
  auto const n_comm_cells = Vector3i{
  dd.cell_grid[1] * dd.cell_grid[2],
  dd.cell_grid[2] * dd.ghost_cell_grid[0],
  dd.ghost_cell_grid[0] * dd.ghost_cell_grid[1]};

  int cnt = 0;
  /* direction loop: x, y, z */
  for (int dir = 0; dir < 3; dir++) {
    lc[(dir + 1) % 3] = 1 - done[(dir + 1) % 3];
    lc[(dir + 2) % 3] = 1 - done[(dir + 2) % 3];
    hc[(dir + 1) % 3] = dd.cell_grid[(dir + 1) % 3] + done[(dir + 1) % 3];
    hc[(dir + 2) % 3] = dd.cell_grid[(dir + 2) % 3] + done[(dir + 2) % 3];
    /* lr loop: left right */
    /* here we could in principle build in a one sided ghost
       communication, simply by taking the lr loop only over one
       value */
    for (int lr = 0; lr < 2; lr++) {
      if (node_grid[dir] == 1) {
        /* just copy cells on a single node */
        if (PERIODIC(dir) || (boundary[2 * dir + lr] == 0)) {
          comm->comm[cnt].type = GHOST_LOCL;
          comm->comm[cnt].node = this_node;

          /* Buffer has to contain Send and Recv cells -> factor 2 */
          comm->comm[cnt].part_lists.resize(2 * n_comm_cells[dir]);
          /* prepare folding of ghost positions */
          if (boundary[2 * dir + lr] != 0) {
              update_component(comm->comm[cnt].shift, boundary[2 * dir + lr] * box_l[dir], dir);
          }

          /* fill send comm cells */
          lc[dir] = hc[dir] = 1 + lr * (dd.cell_grid[dir] - 1);

          dd_fill_comm_cell_lists(comm->comm[cnt].part_lists.data(), lc, hc);
          CELL_TRACE(fprintf(
              stderr,
              "%d: prep_comm %d copy to          grid (%d,%d,%d)-(%d,%d,%d)\n",
              this_node, cnt, lc[0], lc[1], lc[2], hc[0], hc[1], hc[2]));

          /* fill recv comm cells */
          lc[dir] = hc[dir] = 0 + (1 - lr) * (dd.cell_grid[dir] + 1);

          /* place receive cells after send cells */
          dd_fill_comm_cell_lists(
              &comm->comm[cnt].part_lists[n_comm_cells[dir]], lc, hc);
          CELL_TRACE(fprintf(
              stderr,
              "%d: prep_comm %d copy from        grid (%d,%d,%d)-(%d,%d,%d)\n",
              this_node, cnt, lc[0], lc[1], lc[2], hc[0], hc[1], hc[2]));
          cnt++;
        }
      } else {
        /* i: send/recv loop */
        for (i = 0; i < 2; i++) {
          if (PERIODIC(dir) || (boundary[2 * dir + lr] == 0))
            if ((node_pos[dir] + i) % 2 == 0) {
              comm->comm[cnt].type = GHOST_SEND;
              comm->comm[cnt].node = node_neighbors[2 * dir + lr];
              comm->comm[cnt].part_lists.resize(n_comm_cells[dir]);
              /* prepare folding of ghost positions */
              if (boundary[2 * dir + lr] != 0) {
                  update_component(comm->comm[cnt].shift, boundary[2 * dir + lr] * box_l[dir], dir);
              }

              lc[dir] = hc[dir] = 1 + lr * (dd.cell_grid[dir] - 1);

              dd_fill_comm_cell_lists(comm->comm[cnt].part_lists.data(), lc, hc);

              CELL_TRACE(fprintf(stderr,
                                 "%d: prep_comm %d send to   node %d "
                                 "grid (%d,%d,%d)-(%d,%d,%d)\n",
                                 this_node, cnt, comm->comm[cnt].node, lc[0],
                                 lc[1], lc[2], hc[0], hc[1], hc[2]));
              cnt++;
            }
          if (PERIODIC(dir) || (boundary[2 * dir + (1 - lr)] == 0))
            if ((node_pos[dir] + (1 - i)) % 2 == 0) {
              comm->comm[cnt].type = GHOST_RECV;
              comm->comm[cnt].node = node_neighbors[2 * dir + (1 - lr)];
              comm->comm[cnt].part_lists.resize(n_comm_cells[dir]);

              lc[dir] = hc[dir] = (1 - lr) * (dd.cell_grid[dir] + 1);

              dd_fill_comm_cell_lists(comm->comm[cnt].part_lists.data(), lc, hc);
              CELL_TRACE(fprintf(stderr,
                                 "%d: prep_comm %d recv from node %d "
                                 "grid (%d,%d,%d)-(%d,%d,%d)\n",
                                 this_node, cnt, comm->comm[cnt].node, lc[0],
                                 lc[1], lc[2], hc[0], hc[1], hc[2]));
              cnt++;
            }
        }
      }
      done[dir] = 1;
    }
  }
}

/** Revert the order of a communicator: After calling this the
    communicator is working in reverted order with exchanged
    communication types GHOST_SEND <-> GHOST_RECV. */
void dd_revert_comm_order(GhostCommunicator *comm) {
  /* revert order */
  boost::reverse(comm->comm);

  /* exchange SEND/RECV */
  for (auto &c : comm->comm) {
    if(c.shift) {
      c.shift = -c.shift.get();
    }
    if (c.type == GHOST_SEND)
      c.type = GHOST_RECV;
    else if (c.type == GHOST_RECV)
      c.type = GHOST_SEND;
    else if (c.type == GHOST_LOCL) {
      auto first_list = c.part_lists.begin();
      auto last_list = c.part_lists.end();
      auto const n_lists = std::distance(first_list, last_list);
      std::rotate(first_list, first_list + (n_lists / 2), last_list);
    }
  }
}

/** Of every two communication rounds, set the first receivers to prefetch and
 * poststore */
void dd_assign_prefetches(GhostCommunicator *comm) {
    for(auto it = comm->comm.begin(); it != comm->comm.end(); it+=2) {
        if((it->type == GHOST_RECV) and (std::next(it)->type == GHOST_SEND)) {
            it->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
            std::next(it)->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
        }
    }
}

/** update the 'shift' member of those GhostCommunicators, which use
    that value to speed up the folding process of its ghost members
    (see \ref dd_prepare_comm for the original), i.e. all which have
    GHOSTTRANS_POSSHFTD or'd into 'data_parts' upon execution of \ref
    dd_prepare_comm. */
void dd_update_communicators_w_boxl() {
  int cnt = 0;

  /* direction loop: x, y, z */
  for (int dir = 0; dir < 3; dir++) {
    /* lr loop: left right */
    for (int lr = 0; lr < 2; lr++) {
      if (node_grid[dir] == 1) {
        if (PERIODIC(dir) || (boundary[2 * dir + lr] == 0)) {
          /* prepare folding of ghost positions */
          if (boundary[2 * dir + lr] != 0) {
              update_component(cell_structure.local_to_ghost_comm.comm[cnt].shift, boundary[2 * dir + lr] * box_l[dir], dir);
              update_component(cell_structure.ghost_to_local_comm.comm[cnt].shift, -boundary[2 * dir + lr] * box_l[dir], dir);
              update_component(cell_structure.exchange_ghosts_comm.comm[cnt].shift, boundary[2 * dir + lr] * box_l[dir], dir);
          }
          cnt++;
        }
      } else {
        /* i: send/recv loop */
        for (int i = 0; i < 2; i++) {
          if (PERIODIC(dir) || (boundary[2 * dir + lr] == 0))
            if ((node_pos[dir] + i) % 2 == 0) {
              /* prepare folding of ghost positions */
              if (boundary[2 * dir + lr] != 0) {
                  update_component(cell_structure.local_to_ghost_comm.comm[cnt].shift, boundary[2 * dir + lr] * box_l[dir], dir);
                  update_component(cell_structure.ghost_to_local_comm.comm[cnt].shift, -boundary[2 * dir + lr] * box_l[dir], dir);
                  update_component(cell_structure.exchange_ghosts_comm.comm[cnt].shift, boundary[2 * dir + lr] * box_l[dir], dir);
              }
              cnt++;
            }
          if (PERIODIC(dir) || (boundary[2 * dir + (1 - lr)] == 0))
            if ((node_pos[dir] + (1 - i)) % 2 == 0) {
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
void dd_init_cell_interactions() {
  int m, n, o, p, q, r, ind1, ind2;

  for (int i = 0; i < 3; i++) {
    if (dd.fully_connected[i] == true and node_grid[i] != 1) {
      runtimeErrorMsg()
          << "Node grid not compatible with fully_connected property";
    }
  }

  /* loop all local cells */
  DD_LOCAL_CELLS_LOOP(m, n, o) {

    ind1 = get_linear_index(
        m, n, o,
        {dd.ghost_cell_grid[0], dd.ghost_cell_grid[1], dd.ghost_cell_grid[2]});

    std::vector<Cell *> red_neighbors;
    std::vector<Cell *> black_neighbors;

    /* loop all neighbor cells */
    int lower_index[3] = {m - 1, n - 1, o - 1};
    int upper_index[3] = {m + 1, n + 1, o + 1};

    for (int i = 0; i < 3; i++) {
      if (dd.fully_connected[i] == true) {
        lower_index[i] = 0;
        upper_index[i] = dd.ghost_cell_grid[i] - 1;
      }
    }

    for (p = lower_index[2]; p <= upper_index[2]; p++)
      for (q = lower_index[1]; q <= upper_index[1]; q++)
        for (r = lower_index[0]; r <= upper_index[0]; r++) {
          ind2 = get_linear_index(r, q, p,
                                  {dd.ghost_cell_grid[0], dd.ghost_cell_grid[1],
                                   dd.ghost_cell_grid[2]});
          if (ind2 > ind1) {
            red_neighbors.push_back(&cells[ind2]);
          } else {
            black_neighbors.push_back(&cells[ind2]);
          }
        }
    cells[ind1].m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
  }
}

/*************************************************/

/** Returns pointer to the cell which corresponds to the position if
    the position is in the nodes spatial domain otherwise a nullptr
    pointer. */
Cell *dd_save_position_to_cell(const Vector3d &pos) {
  int cpos[3];

  for (int i = 0; i < 3; i++) {
    auto const lpos = pos[i] - my_left[i];
    cpos[i] = static_cast<int>(std::floor(lpos * dd.inv_cell_size[i])) + 1;

    /* particles outside our box. Still take them if
       nonperiodic boundary */
    if (cpos[i] < 1) {
      if (!PERIODIC(i) && boundary[2 * i])
        cpos[i] = 1;
      else
        return nullptr;
    } else if (cpos[i] > dd.cell_grid[i]) {
      if (!PERIODIC(i) && boundary[2 * i + 1])
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

void dd_on_geometry_change(int flags) {
  /* check that the CPU domains are still sufficiently large. */
  for (int i = 0; i < 3; i++)
    if (local_box_l[i] < max_range) {
      runtimeErrorMsg() << "box_l in direction " << i << " is too small";
    }

  /* A full resorting is necessary if the grid has changed. We simply
     don't have anything fast for this case. Probably also not
     necessary. */
  if (flags & CELL_FLAG_GRIDCHANGED) {
    CELL_TRACE(
        fprintf(stderr, "%d: dd_on_geometry_change full redo\n", this_node));

    /* Reset min num cells to default */
    min_num_cells = calc_processor_min_num_cells();

    cells_re_init(CELL_STRUCTURE_CURRENT);
    return;
  }

  /* otherwise, re-set our geometrical dimensions which have changed
     (in addition to the general ones that \ref grid_changed_box_l
     takes care of) */
  for (int i = 0; i < 3; i++) {
    dd.cell_size[i] = local_box_l[i] / (double)dd.cell_grid[i];
    dd.inv_cell_size[i] = 1.0 / dd.cell_size[i];
  }

  double min_cell_size =
      std::min(std::min(dd.cell_size[0], dd.cell_size[1]), dd.cell_size[2]);
  max_skin = min_cell_size - max_cut;

  CELL_TRACE(fprintf(stderr,
                     "%d: dd_on_geometry_change: max_range = %f, "
                     "min_cell_size = %f, max_skin = %f\n",
                     this_node, max_range, min_cell_size, max_skin));

  if (max_range > min_cell_size) {
    /* if new box length leads to too small cells, redo cell structure
       using smaller number of cells. */
    cells_re_init(CELL_STRUCTURE_DOMDEC);
    return;
  }

  /* If we are not in a hurry, check if we can maybe optimize the cell
     system by using smaller cells. */
  if (!(flags & CELL_FLAG_FAST) && max_range > 0) {
    int i;
    for (i = 0; i < 3; i++) {
      int poss_size = (int)floor(local_box_l[i] / max_range);
      if (poss_size > dd.cell_grid[i])
        break;
    }
    if (i < 3) {
      /* new range/box length allow smaller cells, redo cell structure,
         possibly using smaller number of cells. */
      cells_re_init(CELL_STRUCTURE_DOMDEC);
      return;
    }
  }
  dd_update_communicators_w_boxl();
}

/************************************************************/
void dd_topology_init(CellPList *old) {
  int c, p;
  int exchange_data, update_data;

  CELL_TRACE(fprintf(stderr,
                     "%d: dd_topology_init: Number of recieved cells=%d\n",
                     this_node, old->n));

  /* Min num cells can not be smaller than calc_processor_min_num_cells,
     but may be set to a larger value by the user for performance reasons. */
  min_num_cells = std::max(min_num_cells, calc_processor_min_num_cells());

  cell_structure.type = CELL_STRUCTURE_DOMDEC;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = dd_save_position_to_cell;

  /* set up new domain decomposition cell structure */
  dd_create_cell_grid();
  /* mark cells */
  dd_mark_cells();

  /* create communicators */
  dd_prepare_comm(&cell_structure.ghost_cells_comm, GHOSTTRANS_PARTNUM);

  exchange_data =
      (GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  update_data = (GHOSTTRANS_POSITION);

  dd_prepare_comm(&cell_structure.local_to_ghost_comm, 0);
  dd_prepare_comm(&cell_structure.ghost_to_local_comm, 0);
  dd_revert_comm_order(&cell_structure.ghost_to_local_comm);

  dd_assign_prefetches(&cell_structure.local_to_ghost_comm);
  dd_assign_prefetches(&cell_structure.ghost_to_local_comm);

  dd_prepare_comm(&cell_structure.exchange_ghosts_comm, exchange_data);

  dd_assign_prefetches(&cell_structure.ghost_cells_comm);
  dd_assign_prefetches(&cell_structure.exchange_ghosts_comm);

  dd_init_cell_interactions();

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    Particle *part = old->cell[c]->part;
    int np = old->cell[c]->n;
    for (p = 0; p < np; p++) {
      Cell *nc = dd_save_position_to_cell(part[p].r.p);
      /* particle does not belong to this node. Just stow away
         somewhere for the moment */
      if (nc == nullptr)
        nc = local_cells.cell[0];
      append_unindexed_particle(nc, std::move(part[p]));
    }
  }
  for (c = 0; c < local_cells.n; c++) {
    update_local_particles(local_cells.cell[c]);
  }
  CELL_TRACE(fprintf(stderr, "%d: dd_topology_init: done\n", this_node));
}

/************************************************************/
void dd_topology_release() {
}

namespace {
/**
 * @brief Move particles into the cell system if
 *        it belongs to this node.
 *
 * Moves all particles from src into the local cell
 * system if they do belong here. Otherwise the
 * particles are moved into rest.
 *
 * @param src Particles to move.
 * @param rest Output list for left-over particles.
 */
void move_if_local(ParticleList &src, ParticleList &rest) {
  for (int i = 0; i < src.n; i++) {
    auto &part = src.part[i];

    assert(local_particles[src.part[i].p.identity] == nullptr);

    auto target_cell = dd_save_position_to_cell(part.r.p);

    if (target_cell) {
      append_indexed_particle(target_cell, std::move(src.part[i]));
    } else {

      append_unindexed_particle(&rest, std::move(src.part[i]));
    }
  }

  realloc_particlelist(&src, src.n = 0);
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
  for (int i = 0; i < src.n; i++) {
    auto &part = src.part[i];

    assert(local_particles[src.part[i].p.identity] == nullptr);

    if (get_mi_coord(part.r.p[dir], my_left[dir], dir) < 0.0) {
      if (PERIODIC(dir) || (boundary[2 * dir] == 0)) {

        move_unindexed_particle(&left, &src, i);
        if (i < src.n)
          i--;
      }
    } else if (get_mi_coord(part.r.p[dir], my_right[dir], dir) >= 0.0) {
      if (PERIODIC(dir) || (boundary[2 * dir + 1] == 0)) {

        move_unindexed_particle(&right, &src, i);
        if (i < src.n)
          i--;
      }
    }
  }
}

void exchange_neighbors(ParticleList *pl) {
  for (int dir = 0; dir < 3; dir++) {
    /* Single node direction, no action needed. */
    if (node_grid[dir] == 1) {
      continue;
      /* In this (common) case left and right neighbors are
         the same, and we need only one communication */
    } else if (node_grid[dir] == 2) {
      ParticleList send_buf, recv_buf;
      move_left_or_right(*pl, send_buf, send_buf, dir);

      Utils::Mpi::sendrecv(comm_cart, node_neighbors[2 * dir], 0xaa, send_buf,
                           node_neighbors[2 * dir], 0xaa, recv_buf);

      realloc_particlelist(&send_buf, 0);

      move_if_local(recv_buf, *pl);
    } else {
      ParticleList send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;

      move_left_or_right(*pl, send_buf_l, send_buf_r, dir);
      auto reqs_l = Utils::Mpi::isendrecv(
          comm_cart, node_neighbors[2 * dir], 0xaa, send_buf_l,
          node_neighbors[2 * dir], 0xaa, recv_buf_l);
      Utils::Mpi::sendrecv(comm_cart, node_neighbors[2 * dir + 1], 0xaa,
                           send_buf_r, node_neighbors[2 * dir + 1], 0xaa,
                           recv_buf_r);
      boost::mpi::wait_all(reqs_l.begin(), reqs_l.end());

      realloc_particlelist(&send_buf_l, 0);
      realloc_particlelist(&send_buf_r, 0);

      move_if_local(recv_buf_l, *pl);
      move_if_local(recv_buf_r, *pl);
    }
  }
}
} // namespace

void dd_exchange_and_sort_particles(int global, ParticleList *pl) {
  if (global) {
    /* Worst case we need node_grid - 1 rounds per direction.
     * This correctly implies that if there is only one node,
     * no action should be taken. */
    int rounds_left = node_grid[0] + node_grid[1] + node_grid[2] - 3;
    for (; rounds_left > 0; rounds_left--) {
      exchange_neighbors(pl);

      auto left_over =
          boost::mpi::all_reduce(comm_cart, pl->n, std::plus<int>());

      if (left_over == 0) {
        break;
      }
    }
  } else {
    exchange_neighbors(pl);
  }
}

/*************************************************/

int calc_processor_min_num_cells() {
  int i, min = 1;
  /* the minimal number of cells can be lower if there are at least two nodes
     serving a direction,
     since this also ensures that the cell size is at most half the box length.
     However, if there is
     only one processor for a direction, there have to be at least two cells for
     this direction. */
  for (i = 0; i < 3; i++)
    if (node_grid[i] == 1)
      min *= 2;
  return min;
}

/************************************************************/
